#ifndef GEOLOGICAL_MODEL_H
#define GEOLOGICAL_MODEL_H
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

// Represents an isosurface with its associated isovalue and mesh data
struct Isosurface
{
    double isovalue;  // The isovalue defining the isosurface
    std::string mesh; // Placeholder for mesh data (could be a path, ID, or actual mesh object)
};

// Represents a scalar field with its name, values, isosurfaces, and truncation rules
struct ScalarField
{
    std::string name;                                        // Name of the scalar field
    std::vector<double> values;                              // Scalar field values (e.g., sampled grid or function)
    std::vector<Isosurface> isosurfaces;                     // List of isosurfaces for this scalar field
    std::unordered_map<std::string, double> truncationRules; // Truncation rules: other scalar field names and their isovalues
};

// Represents the geological model as a collection of scalar fields and their relationships
class GeologicalModel
{
public:
    // Add a scalar field to the model
    void addScalarField(const ScalarField &field)
    {
        scalarFields.push_back(field);
    }

    // Retrieve a scalar field by name
    ScalarField *getScalarField(const std::string &name)
    {
        for (auto &field : scalarFields)
        {
            if (field.name == name)
            {
                return &field;
            }
        }
        return nullptr; // Not found
    }

    // List all scalar fields
    const std::vector<ScalarField> &getScalarFields() const
    {
        return scalarFields;
    }

private:
    std::vector<ScalarField> scalarFields; // Collection of scalar fields
};

#endif // GEOLOGICAL_MODEL_H