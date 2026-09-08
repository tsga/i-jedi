#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

#include "ijedi_run_add_land_increment.h"

//namespace gdasapp {
  /**
   * AddLandIncrement Class Implementation
   *
   */

  class AddLandIncrement : public oops::Application {
   public:
    explicit AddLandIncrement(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}
    static const std::string classname() {return "AddLandIncrement";}

    int execute(const eckit::Configuration & fullConfig) const {
      // Initialize the application
      // Pass the full configuration to the calculation class
      runAddLandIncrement rinc(fullConfig, this->getComm());
      rinc.run();
      return 0;
    }

   private:
      // -----------------------------------------------------------------------------
      std::string appname() const {
        return "AddLandIncrement";
    }
  };
//}  // namespace gdasapp
