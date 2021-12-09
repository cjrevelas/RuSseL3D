# Changelog

All notable changes to this project will be documented in this file. See [standard-version](https://github.com/conventional-changelog/standard-version) for commit guidelines.

## [3.2.0](https://github.com/mokkapps/changelog-generator-demo/compare/v3.1.0...v3.2.0) (2021-12-09)


### Features

* add compute_phi_end_middle subroutine ([27fcf5b](https://github.com/mokkapps/changelog-generator-demo/commits/27fcf5b7eb6b800fdb7893c24772df70e5e7cee2))
* add exports, computes, gets and binning histogram ([47d2933](https://github.com/mokkapps/changelog-generator-demo/commits/47d2933f7c88d3e941aae758cda9e9ae729587d6))
* add flag enum module ([e45b18a](https://github.com/mokkapps/changelog-generator-demo/commits/e45b18af8cd5d90d4b9bd1dec93b8dca5d2d9463))
* add hybrid contour discr scheme ([5657d2a](https://github.com/mokkapps/changelog-generator-demo/commits/5657d2a381673d1add2d12c56a2b842184a3fe6d))
* add phi_total variable ([127e84b](https://github.com/mokkapps/changelog-generator-demo/commits/127e84bee1747d1a3da000f40953913948bfdb93))
* add square-gradient term (experimental) ([c967711](https://github.com/mokkapps/changelog-generator-demo/commits/c967711bcb46a639a2be55d0166a53589944da66))
* add useful unit transformation constants ([a8ad62c](https://github.com/mokkapps/changelog-generator-demo/commits/a8ad62c8b32f2f2ece682b86ff44c77ef30055ab))
* choose whether to export stuff from input ([259cb05](https://github.com/mokkapps/changelog-generator-demo/commits/259cb053dff7fe6a5e9949f3dce0c40703fbd4a4))
* **constants:** add some useful constants ([3e3d8a5](https://github.com/mokkapps/changelog-generator-demo/commits/3e3d8a59d12d585d9f7bc36633283e3832f2910f))
* contour discr options are: unif, symm, asymm, hybrid ([6e62818](https://github.com/mokkapps/changelog-generator-demo/commits/6e62818d63091059dd2f9f6fcc03ea60d5da4096))
* **energies:** normalize term4 with respect to distance ([cbd41be](https://github.com/mokkapps/changelog-generator-demo/commits/cbd41be72433922b45455bdef60665ae11eb8e13))
* generalize nonbonded interactions and add SL EoS ([b5b1131](https://github.com/mokkapps/changelog-generator-demo/commits/b5b113193cdaedc386ed91d66a3d44dde3556370))
* **main:** add free energy convergence criterion ([abe680a](https://github.com/mokkapps/changelog-generator-demo/commits/abe680a931fbb9be0b1988f4cf185ac38d918595))
* **mesh:** compute number of elements per node ([f9568cd](https://github.com/mokkapps/changelog-generator-demo/commits/f9568cd1c4092dc031646a4540b612e4e4a0f5a7))
* read contour step, dN, instead of number of steps, n ([2629815](https://github.com/mokkapps/changelog-generator-demo/commits/26298151dd8f86e572d79bf97ea022fe1725595a))
* separate io files in debug (d.) and output (o.) ([c9c3f30](https://github.com/mokkapps/changelog-generator-demo/commits/c9c3f30f7149963473d437c6dfb3ad5baf075c03))


### Bug Fixes

* **compute_phi_end_middle:** switch from 1D to 3D ([b30b858](https://github.com/mokkapps/changelog-generator-demo/commits/b30b858431d692ae75100d4be4d6089ac1cb0ae2))
* **eos:** check whether rho_tilde is greater than 1.0 ([414a806](https://github.com/mokkapps/changelog-generator-demo/commits/414a80650d1c7465155029d9086110ac1b212a10))

## [3.1.0](https://github.com/mokkapps/changelog-generator-demo/compare/v3.0.0...v3.1.0) (2021-12-06)


### Features

* add calc_delta_every flag ([366488e](https://github.com/mokkapps/changelog-generator-demo/commits/366488ee3408821e142ef5c41f3dfe752b8eb805))
* add hard sphere wall for Hamaker (wall_dist) ([a63b156](https://github.com/mokkapps/changelog-generator-demo/commits/a63b15691852dad79e26e70c7833dbcb12ae2375))
* contour-step is a function of the current step ([2546c35](https://github.com/mokkapps/changelog-generator-demo/commits/2546c35608521d8a60bb6745fd86735968a8f779))
* determine numerical and analytic delta for grafted chains ([6ce7c02](https://github.com/mokkapps/changelog-generator-demo/commits/6ce7c02b83f8036420d27322548a1e2f4871dc33))
* **energies:** export mx propagator and term4 ([c611c3b](https://github.com/mokkapps/changelog-generator-demo/commits/c611c3b58348818f643b8409fb18a7ba6d03917c))
* **force_fields:** add sphere-plate Hamaker potential ([ca9fd5e](https://github.com/mokkapps/changelog-generator-demo/commits/ca9fd5edd7741dc431145799c670fdd47971310e))
* **iofiles:** gathers all io file names ([b802a56](https://github.com/mokkapps/changelog-generator-demo/commits/b802a56db115691ba2193aa815f91ed706c448ee))
* **mesh:** upgrade mesh parser to be more flexible ([41bdcb9](https://github.com/mokkapps/changelog-generator-demo/commits/41bdcb9ef591aace9a5abde40603a7cd355e2c33))
* set default hard sphere wall dist equal to 5.0 Angstrom ([f7d741a](https://github.com/mokkapps/changelog-generator-demo/commits/f7d741a6efe56cbd3de0513212241330135d541d))


### Bug Fixes

* **adh_ten:** change sign ([cf5823c](https://github.com/mokkapps/changelog-generator-demo/commits/cf5823c1062a114fc43b7027bc889c9773252020))
* **force_fields:** h_12 input argument is not initialized to 0.0 ([a837d29](https://github.com/mokkapps/changelog-generator-demo/commits/a837d2952431c903eca7aa9a1f4d5651a3da19d9))
* **get_sys_time:** remove junk time.out.txt file ([5c5957e](https://github.com/mokkapps/changelog-generator-demo/commits/5c5957e01f1df6008355ad4b594e34d5a4182e47))
* **init_field:** determine hss by removing wall_dist ([cf80383](https://github.com/mokkapps/changelog-generator-demo/commits/cf80383f0856056ac6654cd6ad445186ff88613f))
* **parser:** read input variables without format ([f4a4333](https://github.com/mokkapps/changelog-generator-demo/commits/f4a4333f60233446382c107a7e938e5ea41e1ec1))

## [3.0.0](https://github.com/mokkapps/changelog-generator-demo/compare/v2.2.0...v3.0.0) (2021-12-04)


### Features

* add chain contour interpolation functionality ([d4f7e3a](https://github.com/mokkapps/changelog-generator-demo/commits/d4f7e3ab42acce9611e19784a24821a94c175bc5))
* add functionality for spherical nanoparticles ([a43963a](https://github.com/mokkapps/changelog-generator-demo/commits/a43963ae0077314de06b8a7f2db3bf89bb751e4a))
* add generic force_fields module to replace surf_pot ([7bceb53](https://github.com/mokkapps/changelog-generator-demo/commits/7bceb5399f559b094451c3e386a35a87b7fcd6d4))
* add geometry module to hold geom/mesh variables ([cdcf120](https://github.com/mokkapps/changelog-generator-demo/commits/cdcf1204a7b227af0eee13c6aa3ef6cd4eccca88))
* add option for nonuniform chain contour stepping ([96f295d](https://github.com/mokkapps/changelog-generator-demo/commits/96f295d68ad285c2035351b3a9e06e45838de1fe))
* automatic detection of Dirichlet face ids ([ffb1b88](https://github.com/mokkapps/changelog-generator-demo/commits/ffb1b8893bac2294ec839c5c8b7583ffa9720122))
* determine contour value at all contour points ([212b14f](https://github.com/mokkapps/changelog-generator-demo/commits/212b14f696a5f86b7f3525794a8d336715366fa7))
* different ns for diffusion and convolution ([211ffe9](https://github.com/mokkapps/changelog-generator-demo/commits/211ffe93552d0d210d250cb35bdeac3aa62a144d))
* export all individual energy terms ([224e803](https://github.com/mokkapps/changelog-generator-demo/commits/224e80391a036d3258615d84f7e9d72c8014a55e))
* interf_area is read from parser ([3dc37c9](https://github.com/mokkapps/changelog-generator-demo/commits/3dc37c96f443ccc07d825b4584df1188ec564c31))
* **mesh:** enable auto nanoparticle face ids detection ([ea2c868](https://github.com/mokkapps/changelog-generator-demo/commits/ea2c868d05f6e2cc4632dc601b6493fd8c8a1faa))
* read chain length and ns instead of ds and ns ([b51c007](https://github.com/mokkapps/changelog-generator-demo/commits/b51c007feccfefd202c34e501bf8b88474320390))
* **surf_pot:** export attr, rep and total surface potential ([3bc7925](https://github.com/mokkapps/changelog-generator-demo/commits/3bc7925b7335aa141368fd16418f41fb1841310e))


### Bug Fixes

* **adh_ten:** correct calculation of term4 ([5fde860](https://github.com/mokkapps/changelog-generator-demo/commits/5fde860b3e58885085feaada9ba7681ce4a5e873))
* **adh_ten:** remove the doubling of interf_area ([5b7d456](https://github.com/mokkapps/changelog-generator-demo/commits/5b7d456f0f0442090b42bde016f44fe9c322e662))
* **adh_ten:** term4 takes into account all grafted chains ([df586b8](https://github.com/mokkapps/changelog-generator-demo/commits/df586b8a0e70ad78cde0599a2ccdfcf2942e6b22))
* **grafted_init_cond:** read number of gps without format ([3f4bd43](https://github.com/mokkapps/changelog-generator-demo/commits/3f4bd436599bdc3dadb31bf7c1adadc5e43ce375))
* temporary fix for distance on 3D Hamaker potantial ([d44ea41](https://github.com/mokkapps/changelog-generator-demo/commits/d44ea410f78f6aa220665fbf056be52f40accddf))

## [2.2.0](https://github.com/mokkapps/changelog-generator-demo/compare/v2.1.0...v2.2.0) (2021-12-01)


### Features

* add a subroutine to measure system time ([4fbbd65](https://github.com/mokkapps/changelog-generator-demo/commits/4fbbd6546a4b6d12ba64647842a0c2f331cf6781))
* Edwards is solved with Rg2_per_mon ([4415645](https://github.com/mokkapps/changelog-generator-demo/commits/4415645c58cd13b2930b6d7b9aeec07a858c7473))
* Makefile organizes code in src, obj and run dirs ([c83b753](https://github.com/mokkapps/changelog-generator-demo/commits/c83b753bf822700bd02590e458f09db3a3339d86))
* matrix and grafted chains can have different chain length ([85daf01](https://github.com/mokkapps/changelog-generator-demo/commits/85daf0104ad61506f04032255e73f17779f104b1))
* module files are contained in mod directory ([f5500af](https://github.com/mokkapps/changelog-generator-demo/commits/f5500afff5c40188a0ce62d486e23b76f45758ca))
* option to solve for grafted is read from input ([e546a64](https://github.com/mokkapps/changelog-generator-demo/commits/e546a640f3a412c22811db4d36f26883f9bd848c))
* **qprint:** exports propagators of mx and gr chains ([c3c37e1](https://github.com/mokkapps/changelog-generator-demo/commits/c3c37e14b85d231a2277177d61912891c4f2e63c))
* solve Edwards for grafted chains ([0fe8041](https://github.com/mokkapps/changelog-generator-demo/commits/0fe8041440aea61bb8615760836c181b5d6e8951))
* solve Edwards in presence of multiple grafting points ([00a23d2](https://github.com/mokkapps/changelog-generator-demo/commits/00a23d28a40d8559f7097e62a6486607b96f858f))


### Bug Fixes

* **scfinout:** does not crash when dir faces do not exist ([b6d564a](https://github.com/mokkapps/changelog-generator-demo/commits/b6d564a7d02c522e27bededb755476d386a358c0))
* **tetsh:** correct linear shape functions ([b4ad0ed](https://github.com/mokkapps/changelog-generator-demo/commits/b4ad0ed2d9077362636e7c46edbd71f05423571b))

## [2.1.0](https://github.com/mokkapps/changelog-generator-demo/compare/v2.0.0...v2.1.0) (2021-11-27)


### Features

* add error handling and enhance logging ([da46901](https://github.com/mokkapps/changelog-generator-demo/commits/da46901849069dbfb457889913287ec1acc5e180))
* add fraction variable for field mixing ([bb7ea1e](https://github.com/mokkapps/changelog-generator-demo/commits/bb7ea1e5102be31403a68ef8e228367d5f2d5a1b))
* add write_helper sub ([e3643e3](https://github.com/mokkapps/changelog-generator-demo/commits/e3643e38b7feb617b8256240d268f484eb758f14))
* additional schemes to init field and converge ([c42bdc9](https://github.com/mokkapps/changelog-generator-demo/commits/c42bdc9c16393fe5372a773b78e0294c68fa0a35))
* calculate field error statistics ([972f696](https://github.com/mokkapps/changelog-generator-demo/commits/972f696194209b832ea07db8976a063fba80b89d))
* **edwards:** if symm matrix, remove zeros lines and cols ([9256fb7](https://github.com/mokkapps/changelog-generator-demo/commits/9256fb7d6bca9d7b597ef35414cac1a3f434bb52))
* export field in ascii and binary files ([c450744](https://github.com/mokkapps/changelog-generator-demo/commits/c4507444f2d9b61c150fc218f457750ad292fbe7))
* **main:** periodic profiler ([5e1da26](https://github.com/mokkapps/changelog-generator-demo/commits/5e1da260402a57de4b9f40f1283943cfaef32ffd))
* **mesh_io_3d:** retrieve box dimensions from mesh ([6c12d2c](https://github.com/mokkapps/changelog-generator-demo/commits/6c12d2c573fbca1af672b047891cbbc15592cffc))
* **mesh_io:** export mesh profile along z-axis ([5ddc055](https://github.com/mokkapps/changelog-generator-demo/commits/5ddc055aac0edb5be9aefbe5af278387afed3acf))
* option to deal with symmetric matrices and print A_full ([5c3ef61](https://github.com/mokkapps/changelog-generator-demo/commits/5c3ef61c9c034307c93d1a679c0774e5d3e0dfce))
* option to read field from ascii file ([4ba8c59](https://github.com/mokkapps/changelog-generator-demo/commits/4ba8c59e837fa64a2ce7da0456a28a10587cacf7))
* option to reduce field with chain length ([bdbcee8](https://github.com/mokkapps/changelog-generator-demo/commits/bdbcee825ed732fc8a54eed17e827a3fa8c8d41b))
* option to solve PDE with variable ds ([a5a0a36](https://github.com/mokkapps/changelog-generator-demo/commits/a5a0a3685a3c2d775a73828861eecc807d83d6ce))
* print progress percentage of Edwards solution ([162e762](https://github.com/mokkapps/changelog-generator-demo/commits/162e762c426b090d91da9c10819bff27fe0b8353))
* read ids of dirichlet faces from input ([f101353](https://github.com/mokkapps/changelog-generator-demo/commits/f101353c80686eb955e41c432da99b0757d76681))
* read prof dimension and divide surf ten by area ([24f4e3b](https://github.com/mokkapps/changelog-generator-demo/commits/24f4e3bb938c0c3ec08e1a63c6348b7c7f475f80))
* **scfinout:** input file is now read by parser ([a65b9b7](https://github.com/mokkapps/changelog-generator-demo/commits/a65b9b70670eadfa52047d950d65590bb470342d))
* **spat3d:** warning when mesh and box volume differ ([8ef0363](https://github.com/mokkapps/changelog-generator-demo/commits/8ef0363ae8293f926991dbba164c49f76422508d))


### Bug Fixes

* **adh_ten:** term1, term2 have proper units (m3) ([bdfed83](https://github.com/mokkapps/changelog-generator-demo/commits/bdfed83c6cc8ebdc70406b105c75227a8838ba77))
* last time-step of qf is ns+1 ([028f970](https://github.com/mokkapps/changelog-generator-demo/commits/028f970b13407bc44c062f7602d20dbf4bffd49b))
* **main:** remove integer to real(8) comparison ([324ad68](https://github.com/mokkapps/changelog-generator-demo/commits/324ad68b0910c95aa929f94aa99a06c707ce725c))
* **matrix_assemble:** initialize k, c, w matrices ([d2663c3](https://github.com/mokkapps/changelog-generator-demo/commits/d2663c3e55192840102108745d8d084f6b25c600))
* **MPI:** slave procs receive MUMPS options from master ([f338638](https://github.com/mokkapps/changelog-generator-demo/commits/f338638e5d8c2660e0b220aac03572f302a9a496))
* **spat_3d:** remove u_spat=5. assignment ([0aa0076](https://github.com/mokkapps/changelog-generator-demo/commits/0aa00763115647b61155727303082e773a029195))
* **tetsh:** inequality comparison of real variable ([24622cb](https://github.com/mokkapps/changelog-generator-demo/commits/24622cb2b3e1d97d3594372d219d6b7ea13b0b05))
* **tetshp:** change Jacobian sign from pos to neg ([938e7d5](https://github.com/mokkapps/changelog-generator-demo/commits/938e7d52edb7d25922850c423dfd45bf202a00e6))
* volume is measured in m3 instead of m ([82e5080](https://github.com/mokkapps/changelog-generator-demo/commits/82e5080582c27b502af6e7450521a212a18a6cb8))

## [2.0.0](https://github.com/mokkapps/changelog-generator-demo/compare/v1.0.0...v2.0.0) (2021-11-24)


### Features

* add MPI functionality ([3375f4f](https://github.com/mokkapps/changelog-generator-demo/commits/3375f4f93caf16a1f5c5197cb61d4989f066f242))


### Bug Fixes

* allocate propagator with size ns+1 ([8fc877b](https://github.com/mokkapps/changelog-generator-demo/commits/8fc877bf028e9f8056d0150f894907b59d79ad2a))
* calculate energy terms with 3d integration ([2cba2fb](https://github.com/mokkapps/changelog-generator-demo/commits/2cba2fb4e6d2354237652417fcc686c63b80c312))

## [1.0.0](https://github.com/mokkapps/changelog-generator-demo/compare/v0.0.0...v1.0.0) (2021-11-23)


### Features

* solve 3d edwards equation using 2d matrices ([089a985](https://github.com/mokkapps/changelog-generator-demo/commits/089a98569c5063f78d632c213e8424494656bc0e))

## 0.0.0 (2021-11-18)
