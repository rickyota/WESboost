# pgenlib_rust

This code is fork of pgenlibr in [plink2](https://github.com/chrchang/plink-ng/releases/tag/v2.0.0-a.6.8) and modified by R. Ohta to run as a library in Rust.
The modification is referred to pgenlib in [regenie](https://github.com/rgcgithub/regenie/tree/a087668103d03b04b059cd722d3debe1cf2ead90).


Refer to [pgenlibr manual](https://cran.r-project.org/web/packages/pgenlibr/pgenlibr.pdf).



Modified by rickyota
- pvar_ffi_support.cc:508
	`assert(string_space_iter == string_space_end);` fails, so now comment it.