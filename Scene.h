#pragma once
#ifndef __Scene__
#define __Scene__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

//class Cloth;
class Particle;
class MatrixStack;
class Program;
class Shape;
class RigidBody;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;
	
	double getTime() const { return t; }
	
private:
	double t;
	double h;
	Eigen::Vector3d grav;
	
	//std::shared_ptr<Shape> sphereShape;
	//std::shared_ptr<Cloth> cloth;
	//std::vector< std::shared_ptr<Particle> > spheres;
	std::shared_ptr<RigidBody> rigidbody;
};

#endif
