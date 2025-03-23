#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>
#include <vector>
#include <limits>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;
class Scene;
// -------------------------------------------------
// Global Variables
// -------------------------------------------------
const int Width = 512;
const int Height = 512;
std::vector<float> OutputImage;
// -------------------------------------------------

class Ray {
public:
    vec3 origin;
    vec3 direction;

    Ray(const vec3& o, const vec3& d) : origin(o), direction(normalize(d)) {}
};

class Camera {
public:
    vec3 eye;
    vec3 u, v, w;
    float l, r, b, t, d;
    int nx, ny;

    Camera(const vec3& e, const vec3& u, const vec3& v, const vec3& w, float l, float r, float b, float t, float d, int nx, int ny)
        : eye(e), u(u), v(v), w(w), l(l), r(r), b(b), t(t), d(d), nx(nx), ny(ny) {
    }

    Ray generateRay(int i, int j) const {
        float u_coord = l + (r - l) * (i + 0.5f) / nx;
        float v_coord = b + (t - b) * (j + 0.5f) / ny;
        vec3 dir = normalize(u * u_coord + v * v_coord - w * d);
        return Ray(eye, dir);
    }
};

class Surface {
public:
    virtual ~Surface() = default;
    virtual bool intersect(const Ray& ray, float& t) const = 0;
    virtual vec3 getColor() const = 0;
    virtual vec3 getNormal(const vec3& P) const = 0;
};

class Plane : public Surface {
public:
    float y;
    vec3 color;

    Plane(float y, const vec3& color) : y(y), color(color) {}

    bool intersect(const Ray& ray, float& t) const override {
        if (ray.direction.y == 0) return false;
        t = (y - ray.origin.y) / ray.direction.y;
        return t >= 0;
    }

    vec3 getColor() const override {
        return color;
    }

    vec3 getNormal(const vec3& /*P*/) const override {
        // 평면의 법선은 일정 (여기서는 위쪽 방향)
        return vec3(0.0f, 1.0f, 0.0f);
    }
};

class Sphere : public Surface {
public:
    vec3 center;
    float radius;
    vec3 color;

    Sphere(const vec3& center, float radius, const vec3& color) : center(center), radius(radius), color(color) {}

    bool intersect(const Ray& ray, float& t) const override {
        vec3 oc = ray.origin - center;
        float a = dot(ray.direction, ray.direction);
        float b = 2.0f * dot(oc, ray.direction);
        float c = dot(oc, oc) - radius * radius;
        float discriminant = b * b - 4 * a * c;
        if (discriminant < 0) {
            return false;
        }
        else {
            t = (-b - sqrtf(discriminant)) / (2.0f * a);
            return t >= 0;
        }


    }

    vec3 getColor() const override {
        return color;
    }

    vec3 getNormal(const vec3& P) const override {
        return normalize(P - center);
    }

};

struct Light {
    vec3 position;
    vec3 intensity; 
};

Light sceneLight = { vec3(0.0f, 5.0f, -5.0f), vec3(1.0f, 1.0f, 1.0f) };


vec3 phongShading(const vec3& P, const vec3& N, const vec3& V, const vec3& objectColor) {
    
    float ka = 0.1f;        
    float kd = 0.7f;       
    float ks = 0.2f;        
    float shininess = 16.0f;

    vec3 L = normalize(sceneLight.position - P);
    vec3 ambient = ka * sceneLight.intensity;
    float diff = max(dot(N, L), 0.0f);
    vec3 diffuse = kd * diff * sceneLight.intensity;
    vec3 R = reflect(-L, N);
    float spec = pow(max(dot(V, R), 0.0f), shininess);
    vec3 specular = ks * spec * sceneLight.intensity;
    vec3 color = (ambient + diffuse + specular) * objectColor;
    return clamp(color, 0.0f, 1.0f);
}



struct HitRecord {
    Surface* object;
    float t;
};


class Scene {
public:
    std::vector<Surface*> objects;

    ~Scene() {
        for (auto object : objects) {
            delete object;
        }
    }

    void addObject(Surface* object) {
        objects.push_back(object);
    }


};

bool intersectScene(const Scene& scene, const Ray& ray, HitRecord& hitRec) {
    float t_min = std::numeric_limits<float>::max();
    bool hit = false;
    hitRec.object = nullptr;
    for (const auto& object : scene.objects) {
        float t;
        if (object->intersect(ray, t) && t < t_min) {
            t_min = t;
            hitRec.object = object;
            hit = true;
        }
    }
    hitRec.t = t_min;
    return hit;
}





void render(const Camera& camera, const Scene& scene) {
    OutputImage.clear();
    OutputImage.resize(Width * Height * 3);


    for (int j = Height - 1; j >= 0; --j) {
        for (int i = 0; i < Width; ++i) {
            Ray ray = camera.generateRay(i, j);
            vec3 pixelColor(0.0f); 
            HitRecord hitRec;
            if (intersectScene(scene, ray, hitRec) && hitRec.object) {
                float t = hitRec.t;
                vec3 P = ray.origin + t * ray.direction; 
                vec3 N = hitRec.object->getNormal(P);      
                vec3 V = normalize(camera.eye - P);     
                vec3 objectColor = hitRec.object->getColor(); 
                pixelColor = phongShading(P, N, V, objectColor);
            }
            int index = ((Height - 1 - j) * Width + i) * 3;
            OutputImage[index + 0] = pixelColor.x;
            OutputImage[index + 1] = pixelColor.y;
            OutputImage[index + 2] = pixelColor.z;
        }
    }
}

void resize_callback(GLFWwindow*, int nw, int nh) {
    glViewport(0, 0, nw, nh);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, static_cast<double>(Width), 0.0, static_cast<double>(Height), 1.0, -1.0);
    OutputImage.reserve(Width * Height * 3);
}

int main(int argc, char* argv[]) {
    GLFWwindow* window;
    if (!glfwInit())
        return -1;
    window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);

    if (!window) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return -1;
    }

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glfwSetFramebufferSizeCallback(window, resize_callback);

    vec3 eye = vec3(0.0f, 0.0f, 0.0f);
    vec3 u = vec3(1.0f, 0.0f, 0.0f);
    vec3 v = vec3(0.0f, 1.0f, 0.0f);
    vec3 w = vec3(0.0f, 0.0f, 1.0f);
    float l = -0.1f, r = 0.1f, b = -0.1f, t = 0.1f, d = 0.1f;
    int nx = 512, ny = 512;
    Camera camera(eye, u, v, w, l, r, b, t, d, nx, ny);

    Scene scene;
    //scene.addObject(new Plane(-2.0f, vec3(0.0f, 1.0f, 0.0f))); // Green plane
    scene.addObject(new Sphere(vec3(-4.0f, 0.0f, -7.0f), 1.0f, vec3(1.0f, 0.0f, 0.0f))); // Red sphere
    scene.addObject(new Sphere(vec3(0.0f, 0.0f, -7.0f), 2.0f, vec3(0.0f, 0.0f, 1.0f))); // Blue sphere
    scene.addObject(new Sphere(vec3(4.0f, 0.0f, -7.0f), 1.0f, vec3(1.0f, 1.0f, 0.0f))); // Yellow sphere

    render(camera, scene);

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
        glfwSwapBuffers(window);
        glfwPollEvents();
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}