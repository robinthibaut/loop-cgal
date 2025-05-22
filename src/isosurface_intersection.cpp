#include "isosurface_intersection.h"
#include "marching_cubes.h"
#include <CGAL/Polygon_mesh_processing/corefinement.h>
Mesh compute_optimized_intersection(
    ScalarFields &scalar_field1,
    ScalarFields &scalar_field2,
    double iso_value1,
    double iso_value2,
    double origin_x, double origin_y, double origin_z,
    double cell_size_x, double cell_size_y, double cell_size_z)
{

    // Extract first isosurface with active cells
    TruncationRules truncation_rules1;
    TruncationRules truncation_rules2;
    MarchingCubes mc1(scalar_field1, iso_value1, truncation_rules1,(origin_x, origin_y, origin_z), cell_size_x);
    Mesh mesh1;
    active_cells_set active_cells1;
    std::tie(mesh1, active_cells1) = mc1.generate_mesh();
    MarchingCubes mc2(scalar_field2, iso_value2, truncation_rules2, Point(origin_x, origin_y, origin_z), cell_size_x);
    Mesh mesh2;
    active_cells_set active_cells2;
    std::tie(mesh2, active_cells2) = mc2.generate_mesh();
    // Extract first isosurface with active cells
    
    active_cells_set common_cells;
    for (const auto &cell : active_cells1)
    {
        if (active_cells2.find(cell) != active_cells2.end())
        {
            common_cells.insert(cell);
        }
    }

    // Compute intersection only in common cells
    Mesh result;

    if (common_cells.empty())
    {
        std::cout << "No intersections found between isosurfaces." << std::endl;
        return result;
    }

    // Extract submeshes for each common cell and compute local intersections
    for (const auto &cell : common_cells)
    {
        // Define cell bounding box
        Point min_corner(origin_x + cell.i * cell_size_x,
                         origin_y + cell.j * cell_size_y,
                         origin_z + cell.k * cell_size_z);
        Point max_corner(origin_x + (cell.i + 1) * cell_size_x,
                         origin_y + (cell.j + 1) * cell_size_y,
                         origin_z + (cell.k + 1) * cell_size_z);

        // Extract submeshes contained in this cell
        Mesh submesh1 = extract_submesh(mesh1, min_corner, max_corner);
        Mesh submesh2 = extract_submesh(mesh2, min_corner, max_corner);

        // Compute local intersection
        Mesh local_intersection;
        bool valid = CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(
            submesh1, submesh2, local_intersection);

        if (valid && !local_intersection.is_empty())
        {
            // Merge local intersection with global result
            merge_meshes(result, local_intersection,result);
        }
    }

    return result;
}
// bool facet_is_within_bounds(
//     const Mesh::Facet &facet,
//     const Point &min_corner,
//     const Point &max_corner)
// {
//     for (auto vertex : facet.vertices()) {
//         const Point& p = vertex->point();
//         if (p.x() < min_corner.x() || p.x() > max_corner.x() ||
//             p.y() < min_corner.y() || p.y() > max_corner.y() ||
//             p.z() < min_corner.z() || p.z() > max_corner.z()) {
//             return false;
//         }
//     }
//     return true;
// }

// Mesh extract_submesh(
//     const Mesh &mesh,
//     const Point &min_corner,
//     const Point &max_corner)
// {
//     Mesh submesh;
//     for (const auto &facet : mesh.faces())
//     {
//         // Check if the facet is within the bounding box
//         if (facet_is_within_bounds(facet, min_corner, max_corner))
//         {
//             submesh.add_face(facet);
//         }
//     }
//     return submesh;
// }