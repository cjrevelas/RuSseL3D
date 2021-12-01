# Changelog

All notable changes to this project will be documented in this file. See [standard-version](https://github.com/conventional-changelog/standard-version) for commit guidelines.

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
