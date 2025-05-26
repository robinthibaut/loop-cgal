#ifndef SURFACE_H
#define SURFACE_H
#include <vector>

class Surface
{
    public:
        Surface() = default;
        Surface(const std::vector<std::array<double, 3>> &vertices,
                const std::vector<std::array<int, 3>> &triangles, std::vector<int> material)
            : vertices_(vertices), triangles_(triangles), material_(material) {}

        const std::vector<std::array<double, 3>> &get_vertices() const { return vertices_; }
        const std::vector<std::array<int, 3>> &get_triangles() const { return triangles_; }
        void set_vertices(const std::vector<std::array<double, 3>> &v) { vertices_ = v; }
        void set_triangles(const std::vector<std::array<int, 3>> &t) { triangles_ = t; }
        void split_material(const std::vector<int> &mask);
    private:
        
        std::vector<std::array<double, 3>> vertices_;
        std::vector<std::array<int, 3>> triangles_;
        std::vector<int> material_;
        std::vector<int> material_ids_;
};
#endif // SURFACE_H