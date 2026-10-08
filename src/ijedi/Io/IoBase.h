#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/OptionalPolymorphicParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

namespace ijedi
{
  class Geometry;

  // -------------------------------------------------------------------------------------------------

  class IoBase : public util::Printable, private boost::noncopyable
  {
   public:
    explicit IoBase(const Geometry &, const eckit::LocalConfiguration &);
    virtual ~IoBase() {}
    void readBase(atlas::FieldSet &) const;
    void writeBase(const atlas::FieldSet &) const;

   private:
    // Child read/write methods
    virtual void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
                      const eckit::LocalConfiguration &) const = 0;
    virtual void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
                       const eckit::LocalConfiguration &) const = 0;

    // Child print method
    virtual void print(std::ostream &) const = 0;

    // Configuration holding the field names as used in the files
    eckit::LocalConfiguration fieldIoNames_;
    // Configuration holding scaling factors used to transition units between forecast model
    // and JEDI
    eckit::LocalConfiguration fieldIoScaling_;
  };

  // -------------------------------------------------------------------------------------------------

  class IoParametersBase : public oops::Parameters
  {
    OOPS_ABSTRACT_PARAMETERS(IoParametersBase, Parameters)
   public:
    oops::OptionalParameter<std::string> filetype{"filetype", this};
    oops::OptionalParameter<eckit::LocalConfiguration> fieldIoNames{"field io names", this};
  };

  // -------------------------------------------------------------------------------------------------

  class IoFactory;

  // -------------------------------------------------------------------------------------------------

  // Replace %{member}% in the io strings with the top-level member that OOPS sets
  eckit::LocalConfiguration swapIoMember(const eckit::Configuration &);

  // -------------------------------------------------------------------------------------------------

  class IoParametersWrapper : public oops::Parameters
  {
    OOPS_CONCRETE_PARAMETERS(IoParametersWrapper, Parameters)
   public:
    oops::OptionalPolymorphicParameter<IoParametersBase, IoFactory> ioParameters{"filetype", this};
  };

  // -------------------------------------------------------------------------------------------------

  class IoFactory
  {
   public:
    static IoBase *create(const Geometry &, const IoParametersBase &params);

    static std::unique_ptr<IoParametersBase> createParameters(const std::string &name);

    static std::vector<std::string> getMakerNames()
    {
      return oops::keys(getMakers());
    }

    virtual ~IoFactory() = default;

   protected:
    explicit IoFactory(const std::string &name);

   private:
    virtual IoBase *make(const Geometry &, const IoParametersBase &) = 0;

    virtual std::unique_ptr<IoParametersBase> makeParameters() const = 0;

    static std::map<std::string, IoFactory *> &getMakers()
    {
      static std::map<std::string, IoFactory *> makers_;
      return makers_;
    }
  };

  // -------------------------------------------------------------------------------------------------

  template <class T>
  class IoMaker : public IoFactory
  {
    typedef typename T::Parameters_ Parameters_;

    IoBase *make(const Geometry &geom, const IoParametersBase &params) override
    {
      return new T(geom, dynamic_cast<const Parameters_ &>(params));
    }

    std::unique_ptr<IoParametersBase> makeParameters() const override
    {
      return std::make_unique<Parameters_>();
    }

   public:
    explicit IoMaker(const std::string &name) : IoFactory(name) {}
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi
