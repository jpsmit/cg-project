class CustomCompZ {
    Mesh mesh;

public:

    CustomCompZ(Mesh& meshh) {
        mesh = meshh;
    }

    bool operator() (Triangle left, Triangle right) const {
        float leftZ = (mesh.vertices[left[0]].p.z + mesh.vertices[left[1]].p.z + mesh.vertices[left[2]].p.z) / 3;
        float rightZ = (mesh.vertices[right[0]].p.z + mesh.vertices[right[1]].p.z + mesh.vertices[right[2]].p.z) / 3;

        return leftZ < rightZ; 
    }
};