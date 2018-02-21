#include <iostream>

#include "Scene.h"
#include "Particle.h"
//#include "Cloth.h"
#include "Shape.h"
#include "Program.h"
#include "RigidBody.h"

using namespace std;
using namespace Eigen;

Scene::Scene() :
	t(0.0),
	h(1e-3),
	grav(0.0, 0.0, 0.0)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 1e-3;
	
	grav << 0.0, -9.8, 0.0;
	
	rigidbody = make_shared<RigidBody>();

	//sphereShape = make_shared<Shape>();
	//sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
	
	//auto sphere = make_shared<Particle>(sphereShape);
	//spheres.push_back(sphere);
	//sphere->r = 0.1;
	//sphere->x = Vector3d(0.0, 0.2, 0.0);
}

void Scene::init()
{
	rigidbody->init();
}

void Scene::tare()
{
	//cloth->tare();
}

void Scene::reset()
{
	t = 0.0;

	//cloth->reset();

}

void Scene::step()
{
	t += h;

	rigidbody->step(h);
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
	
	rigidbody->draw(MV, prog);
}
