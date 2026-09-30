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
| `yaluk.ini` | Main YALUK settings (see below). |
| `yaluk_status.ini` | Output and option flags (see below). |
| `CaseFiles/` | Line geometry (`linea_NNN.txt`), lightning current (`corr_NNNNN.txt`) and miscellaneous parameters (`Miscelaneo.txt`). |
| `*.csv` | Lightning database (strike parameters for each event). |
| `*.pch` | ATP punch files, such as surge-arrester models. |

To run a case, execute `tpbig` from inside the case folder (ATP's `startup` file must be there too):

```bash
cd examples/SingleLine
/path/to/tpbig BOTH test.atp s -r
```

YALUK reads both `.ini` files from the directory where `tpbig` runs. Only the first value on each line is read; the text after it (usually starting with `%`) is a comment, so keep the lines in this order.

### yaluk.ini

```text
CaseFiles % Folder with Config_files
1         % Number of lightning to simulate
5         % Maximum number of lines
1         % Maximum number of conductors
```

| Line | Meaning |
| ---- | ------- |
| 1 | Folder, relative to the case directory, that holds the line and lightning files. If it does not exist, YALUK creates it and asks you to put the files there. |
| 2 | Lightning case number. It selects the current file `corr_NNNNN.txt` (zero-padded to five digits, so `1` reads `corr_00001.txt`). |
| 3 | Number of line sections. YALUK reads `linea_001.txt` up to `linea_NNN.txt`, so it must match the number of YALUK instances in the `.atp` file. |
| 4 | Maximum number of conductors per line. Values above 8 are limited to 8, and 0 is replaced by 4 with a warning. |

If `yaluk.ini` is missing, YALUK prints the path it tried and asks for the name of the configuration file on the console.

### yaluk_status.ini

```text
.False.  % Imprimir_campo  Print the electromagnetic field
.False.  % Imprimir_inf    Print process information
.False.  % Read_campo      Read the electromagnetic field from an external file
.False.  % status_perfil   Use a non-constant line profile
```

Each line is a Fortran logical (`.True.` or `.False.`):

| Flag | When `.True.` |
| ---- | ------------- |
| `Imprimir_campo` | Writes the computed electromagnetic field of each line section to `campoNNN.txt`. |
| `Imprimir_inf` | Prints extra information about the process to the console. |
| `Read_campo` | Reserved for reading the electromagnetic field from an external file. The current code reads this flag but does not use it yet. |
| `status_perfil` | Uses a line profile that is not constant along the line (varying conductor heights). |

The last two lines are optional. If `yaluk_status.ini` is missing, all flags are `.False.`.

## License

YALUK is released under the GNU General Public License v3. See [LICENSE](LICENSE). ATP-EMTP and DISLIN are licensed separately by their owners.
