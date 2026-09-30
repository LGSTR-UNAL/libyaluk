# libyaluk

**YALUK** is a library for **ATP-EMTP** that simulates **lightning-induced overvoltages on overhead distribution lines**.

It is linked into ATP as a *foreign model* (`YALUK_DLL_MODEL`). When the ATP solver runs, YALUK computes the electromagnetic field radiated by a nearby lightning strike, couples it to the line conductors, and returns the induced voltages and currents to the ATP network. This lets you study indirect strikes on realistic distribution circuits with transformers, surge arresters, and other ATP components.

## Repository layout

```text
yaluk/
  *.f                 YALUK Fortran sources (yaluk_init, yaluk_exec, EM field and helper routines)
  makefile            Builds libyaluk.a from the YALUK sources
  makebig             Builds tpbig: compiles the ATP user sources and links them with libyaluk.a
examples/             Sample ATP cases (.atp), lightning databases (.csv) and YALUK config files
LICENSE               GNU GPL v3
```

## Requirements

This repository does **not** include ATP or DISLIN. You must get both yourself:

1. **ATP-EMTP**: ATP is licensed software and is distributed only to registered users. Request access from <https://www.atp-emtp.org/> and download the Linux (gfortran) distribution. You need its user-modifiable sources (`dimdef.f`, `newmods.f`, `comtac.f`, `fgnmod.f`, `usrfun.f`, …) and the precompiled `tpbig.a` library. `tpbig.a` holds the ATP core and is **not** built by this repository: it must come with the ATP distribution.
2. **DISLIN**: ATP's Linux build links against the DISLIN plotting library. Get the build that matches your compiler and architecture from <https://www.dislin.de/> and provide it as `dislin.a` (copy or symlink it under that name if the file is called something like `dislin-11.x.a`).
3. **Toolchain**: `gfortran`, `gcc`, and `make` (OpenMP support is needed for `-fopenmp`).

## Building ATP with YALUK

The build expects the ATP files and DISLIN in one directory (`ATPDIR`). By default this is `../libatp` relative to the `yaluk/` folder:

```text
libatp/            <- ATPDIR: ATP sources (with the edited fgnmod.f), tpbig.a, dislin.a
libyaluk/
  yaluk/           <- makefile, makebig, YALUK sources
```

1. Copy the ATP distribution files and `dislin.a` into `ATPDIR`. Check that `ATPDIR/tpbig.a` and `ATPDIR/dislin.a` exist; otherwise make stops with `No rule to make target '.../tpbig.a'`.
2. Edit `ATPDIR/fgnmod.f` to register YALUK, as described in [Registering YALUK in fgnmod.f](#registering-yaluk-in-fgnmodf).
3. Build ATP with YALUK from the repository root:

   ```bash
   make -f yaluk/makebig
   # or, if ATP lives somewhere else:
   make -f yaluk/makebig ATPDIR=/path/to/atp
   ```

   `makebig` first runs `yaluk/makefile`, which compiles the YALUK sources into `yaluk/libyaluk.a` (only the changed files are recompiled). It then compiles the ATP user sources from `ATPDIR` and links everything with `libyaluk.a`, `tpbig.a` and `dislin.a`. The `tpbig` executable is written to `yaluk/`; the intermediate `.o` files are deleted after linking.

   To build only the YALUK library, run `make -C yaluk`. To remove the build products, run `make -f yaluk/makebig clean` (removes `tpbig` and the ATP objects) and `make -C yaluk clean` (removes `libyaluk.a` and `*.mod`).

If your cases need larger ATP tables, resize ATP with its `vardim` utility before linking, as described in the ATP documentation.

### Registering YALUK in fgnmod.f

`fgnmod.f` is part of the ATP distribution, so it is not included here. ATP uses its `FGNMOD` subroutine to map foreign model names to Fortran routines. Make these two changes in the `FGNMOD` subroutine of `ATPDIR/fgnmod.f` (leave `fgnfun` and the sample routines as they are). Keep a copy of the edited file, because a new ATP version will overwrite it.

1. Register the model name in slot 4 of `refnam`:

   ```fortran
         DATA refnam(4) / 'YALUK_DLL_MODEL' /
   ```

   The name must be uppercase and match the name used in `MODEL ... FOREIGN YALUK_DLL_MODEL` in the `.atp` files.

2. Call the YALUK routines in the `iname.EQ.4` branch:

   ```fortran
         ELSE IF ( iname.EQ.4 ) THEN
          IF (iniflg.EQ.1) THEN
           CALL yaluk_init(xdata, xin, xout, xvar)
          ELSE
           CALL yaluk_exec(xdata, xin, xout, xvar)
          ENDIF
   ```

   `yaluk_init` runs once when the model is initialized and `yaluk_exec` runs at every time step. Both are provided by `libyaluk.a`.

If slot 4 is already used by another foreign model, use any free slot `n` and put the calls in the matching `iname.EQ.n` branch. If no slot is free, increase `namcnt` in `FGNMOD`. Remember that `fgnmod.f` is fixed-form Fortran, so statements must start at column 7 or later.

## Running a case

Each case folder in `examples/` contains:

| File | Purpose |
| ---- | ------- |
| `*.atp` | ATP input deck. It declares `MODEL ... FOREIGN YALUK_DLL_MODEL` and one YALUK instance per line section. |
| `yaluk.ini` | Simulation settings: config folder, number of lightning events, maximum lines and conductors. |
| `yaluk_status.ini` | Flags for printing the EM field and process info, reading the field from an external file, and using a non-constant line profile. |
| `*.csv` | Lightning database (strike parameters for each event). |
| `*.pch` | ATP punch files, such as surge-arrester models. |

To run a case, execute `tpbig` from inside the case folder (ATP's `startup` file must be there too):

```bash
cd examples/SingleLine
/path/to/tpbig BOTH test.atp s -r
```

The `run.sh` scripts in the example folders show the full sequence (build, `vardim`, copy `startup`, run), but they use the paths of an older Vagrant setup (`/vagrant/...`); adjust them to your directories before using them.

## License

YALUK is released under the GNU General Public License v3. See [LICENSE](LICENSE). ATP-EMTP and DISLIN are licensed separately by their owners.
