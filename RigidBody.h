#pragma once
#ifndef __RigidBody__
#define __RigidBody__

#include <vector>
#include <memory>
#include <tetgen.h>
#include "QuadProgMosek.h"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class MatrixStack;
class Program;
class Particle;

typedef Eigen::Triplet<double> ETriplet;

class RigidBody {
public:
	RigidBody();
	void init();
	void updatePosNor(Eigen::MatrixXd E);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> p)const;
	Eigen::Vector3d local2world(Eigen::MatrixXd E, Eigen::Vector3d x); // compute the world position given a local position x on a rigid body
	Eigen::Matrix3d vec2crossmatrix(Eigen::Vector3d a); // repackage a vector into a cross-product matrix
	void initParam();
	void computeE();
	void computeB();
	void computeEE();
	void computePHI();
	void computePhiT();
	void computeOMEGA();
	void updateLocalFrame();
	void step(double h);
	tetgenio in, out;
	
	int nVerts;
	int nTriFaces;
	int nEdges;
	double mass;
	int numCol;
	std::vector <std::shared_ptr<Particle>>nodes;

	Eigen::MatrixXd M; // fixed
	

	Eigen::VectorXd Phi; // 6x1 angular and linear velocitie repackeaged into 6-vector
	Eigen::MatrixXd PHI; //
	Eigen::MatrixXd PhiT; //6x6 spatial cross product matrix cosisting of [omega_i]
	Eigen::VectorXd B; // 6x1 forces
	Eigen::MatrixXd E; // 4x4 transformation matrix consisting of rotational and translational components
	Eigen::MatrixXd EE; // 4x4 spatial velocity matrix
	Eigen::MatrixXd Etemp;

	Eigen::Vector3d Omega;
	Eigen::Matrix3d OMEGA; // 3X3 the crossproduct matrix of Omega
	Eigen::Matrix3d VC; // 3x3 the crossproduct matrix of V
	Eigen::Vector3d V;
	Eigen::VectorXd v;
	Eigen::VectorXd f;

	Eigen::Vector3d p; // 3x1 position of the local frame's origin in world coordinates
	Eigen::Matrix3d R; // 3x3 rotation matrix of which each colum corresponds to the frame's basis vectors e_k, expressed in world coordinates
	
	Eigen::Vector3d ynormal;
	Eigen::MatrixXd temp; // 3x6 

	std::shared_ptr<QuadProg> program;
	
	std::vector<ETriplet> A_;
	Eigen::SparseMatrix<double> A;
	std::vector<ETriplet> C_;
	std::vector<int> index;
	Eigen::SparseMatrix<double> C;
	Eigen::VectorXd xl;
	Eigen::VectorXd xu;
	Eigen::VectorXd b;
	Eigen::Matrix3d I;
	Eigen::Vector3d g; // gravity
	Eigen::VectorXd RHS;
	Eigen::VectorXd sol;

private:
	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;
	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;

};

#endif