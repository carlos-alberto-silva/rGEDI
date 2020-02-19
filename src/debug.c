void initR() {
 char *argv[] = {"REmbeddedPostgres", "--gui=none", "--silent"};
 int argc = sizeof(argv)/sizeof(argv[0]);

  Rf_initEmbeddedR(argc, argv);
}