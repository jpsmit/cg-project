
class CustomCompY {
    Mesh mesh;

public:

    CustomCompY(Mesh& meshh) {
        mesh = meshh;
    }

    bool operator() (Triangle left, Triangle right) const {
        //float leftX = std::min(mesh.vertices[left[0]].p.x, std::min(mesh.vertices[left[1]].p.x, mesh.vertices[left[2]].p.x));
        //float rightX = std::min(mesh.vertices[right[0]].p.x, std::min(mesh.vertices[right[1]].p.x, mesh.vertices[right[2]].p.x));

        float leftY = (mesh.vertices[left[0]].p.y + mesh.vertices[left[1]].p.y + mesh.vertices[left[2]].p.y) / 3;
        float rightY = (mesh.vertices[right[0]].p.y + mesh.vertices[right[1]].p.y + mesh.vertices[right[2]].p.y) / 3;

        return leftY < rightY;
    }
};