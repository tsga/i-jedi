#include <unordered_set>

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/IoBase.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  IoFactory::IoFactory(const std::string &name)
  {
    if (getMakers().find(name) != getMakers().end())
    {
      oops::Log::error() << name << " already registered in ijedi::IoFactory." << std::endl;
      ABORT("Element already registered in ijedi::IoFactory.");
    }
    getMakers()[name] = this;
  }

  // -------------------------------------------------------------------------------------------------

  IoBase *IoFactory::create(const Geometry &geom, const IoParametersBase &params)
  {
    oops::Log::trace() << "IoBase::create starting" << std::endl;
    const std::string &id = params.filetype.value().value();
    typename std::map<std::string, IoFactory *>::iterator jloc = getMakers().find(id);
    if (jloc == getMakers().end())
    {
      oops::Log::error() << id << " does not exist in ijedi::IoFactory." << std::endl;
      ABORT("Element does not exist in ijedi::IoFactory.");
    }
    IoBase *ptr = jloc->second->make(geom, params);
    oops::Log::trace() << "IoBase::create done" << std::endl;
    return ptr;
  }

  // -------------------------------------------------------------------------------------------------

  std::unique_ptr<IoParametersBase>
  IoFactory::createParameters(const std::string &name)
  {
    typename std::map<std::string, IoFactory *>::iterator it =
        getMakers().find(name);
    if (it == getMakers().end())
    {
      throw std::runtime_error(name + " does not exist in ijedi::IoFactory");
    }
    return it->second->makeParameters();
  }

  // -------------------------------------------------------------------------------------------------

  IoBase::IoBase(const Geometry &geom, const eckit::LocalConfiguration &conf)
  {
    oops::Log::trace() << "IoBase::IoBase starting" << std::endl;
    // If conf has 'field io names' then extract from the config and overwrite fieldIoNames_
    if (conf.has("field io names"))
    {
      fieldIoNames_ = conf.getSubConfiguration("field io names");
    }
    // If conf has 'field io scaling' then extract from the config and overwrite fieldIoScaling_
    if (conf.has("field io scaling"))
    {
      fieldIoScaling_ = conf.getSubConfiguration("field io scaling");
    }

    // Check the configs for long name correctness
    const std::vector<std::string> fmdLongNames = geom.getFieldMetadata().getLongNames();
    std::unordered_set<std::string> validNames(fmdLongNames.begin(), fmdLongNames.end());

    // Get the keys from the configs
    const std::vector<std::string> fieldIoNamesKeys = fieldIoNames_.keys();
    const std::vector<std::string> fieldIoScalingKeys = fieldIoScaling_.keys();

    // Check for key validity
    for (const std::string &key : fieldIoNamesKeys)
    {
      const std::string msg = "The \"field io names\" configuration contains \"" + key +
                              "\", which is not part of the field metadata.";
      ASSERT_MSG(validNames.find(key) != validNames.end(), msg);
    }
    for (const std::string &key : fieldIoScalingKeys)
    {
      const std::string msg = "The \"field io scaling\" configuration contains \"" + key +
                              "\", which is not part of the field metadata.";
      ASSERT_MSG(validNames.find(key) != validNames.end(), msg);
    }

    oops::Log::trace() << "IoBase::IoBase done" << std::endl;
  }

  // -------------------------------------------------------------------------------------------------

  void IoBase::readBase(atlas::FieldSet &x) const
  {
    // Call read method from the child class
    this->read(x, fieldIoNames_, fieldIoScaling_);
  }

  // -------------------------------------------------------------------------------------------------

  void IoBase::writeBase(const atlas::FieldSet &x) const
  {
    // Call write method from the child class
    this->write(x, fieldIoNames_, fieldIoScaling_);
  }

  // -------------------------------------------------------------------------------------------------

  void IoBase::writeBaseWithConfig(const atlas::FieldSet &x, const eckit::LocalConfiguration &config) const
  {
    // Call write method from the child class
    this->writeWithConfig(x, fieldIoNames_, fieldIoScaling_, config);
  }

  // -------------------------------------------------------------------------------------------------
 
  // Default implementation — just calls write()
  void IoBase::writeWithConfig(const atlas::FieldSet & x,
                             const eckit::LocalConfiguration & fileionames,
                             const eckit::LocalConfiguration & fileioscaling,
                             const eckit::LocalConfiguration & config) const {
    this->write(x, fileionames, fileioscaling);
  }

}  // namespace ijedi
