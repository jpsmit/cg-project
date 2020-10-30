
class CustomCompX {
    Mesh mesh;

public:
    
    CustomCompX(Mesh& meshh) {
        mesh = meshh;
    }

    bool operator() (Triangle left, Triangle right) const {
        float leftX = (mesh.vertices[left[0]].p.x + mesh.vertices[left[1]].p.x + mesh.vertices[left[2]].p.x)/3;
        float rightX = (mesh.vertices[right[0]].p.x + mesh.vertices[right[1]].p.x + mesh.vertices[right[2]].p.x)/3;

        return leftX < rightX;  //mesh.vertices[left[0]].p.z < mesh.vertices[right[0]].p.z;      //leftX < rightX;
    }
};