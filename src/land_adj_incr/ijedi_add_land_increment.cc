#include "ijedi_add_land_increment.h"

#include "ijedi/Utilities/Traits.h"

#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  //AddLandIncrement<ijedi::Traits> addLandIncrement;
  AddLandIncrement addLandIncrement;
  return run.execute(addLandIncrement);
}
