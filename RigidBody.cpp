#define TETLIBRARY
#include "RigidBody.h"
#include <iostream>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Particle.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"
#include <unsupported/Eigen/src/MatrixFunctions/MatrixExponential.h>

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> ETriplet;
double inf = numeric_limits<double>::infinity();

RigidBody::RigidBody() {
	
	ynormal << 0.0, 1.0, 0.0;
	g << 0, -9.8, 0;
	xl.resize(6);
	xu.resize(6);
	xl << -inf, -inf, -inf, -inf, -inf, -inf;
	xu << inf, inf, inf, inf, inf, inf;

	I.setIdentity();

	in.load_ply("cube");
	tetrahedralize("pqz", &in, &out);
	nVerts = out.numberofpoints;
	nTriFaces = out.numberoftrifaces;
	nEdges = out.numberofedges;

	mass = 0.1 * nVerts;

	for (int i = 0; i < nVerts; i++) {
		auto p = make_shared<Particle>();
		nodes.push_back(p);
		// positon in the local frame
		p->x0 << out.pointlist[3 * i], out.pointlist[3 * i + 1], out.pointlist[3 * i + 2];
		p->v << 0.0, 0.0, 0.0;
		p->m = 0.1;
		p->i = i;
	}
	
	// Build vertex buffers
	posBuf.clear();
	norBuf.clear();
	texBuf.clear();
	eleBuf.clear();

	posBuf.resize(nTriFaces * 9);
	norBuf.resize(nTriFaces * 9);
	eleBuf.resize(nTriFaces * 3);

	M.resize(6, 6);
	M.setZero();

	E.resize(4, 4);
	E.setZero();

	EE.resize(4, 4);
	EE.setZero();

	Phi.resize(6);
	Phi.setZero();

	PHI.resize(6, 6);
	PHI.setZero();

	PhiT.resize(6, 6);
	PhiT.setZero();

	B.resize(6);
	B.setZero();

	v.resize(nVerts * 3);
	v.setZero();
	f.resize(nVerts * 3);
	f.setZero();

	temp.resize(3, 6);
	temp.setZero();
	temp.block(0, 3, 3, 3) = I;
	sol.resize(6);
	initParam();

	updatePosNor(E);

	for (int i = 0; i < nTriFaces; i++) {
		for (int j = 0; j < 3; j++) {
			eleBuf[3 * i + j] = 3 * i + j;
		}
	}
}

void RigidBody::initParam() {

	this->M << mass/6.0, 0, 0, 0, 0, 0,
				0, mass/6.0, 0, 0, 0, 0,
				0, 0, mass/6.0, 0, 0, 0,
				0, 0, 0, mass, 0, 0,
				0, 0, 0, 0, mass, 0,
				0, 0, 0, 0, 0, mass;

	for (int t = 0; t < 3; t++) {
		A_.push_back(ETriplet(t, t, mass / 6.0));
		A_.push_back(ETriplet(t+3, t+3, mass));
	}
	
	A.resize(6, 6);
	A.setFromTriplets(A_.begin(), A_.end());

	// set the rotation matrix of the local fram
	this->R = I; // no rotation
	double thetaz = 30.0;
	double thetax = 40.0;
	double thetay = 20.0;
	Matrix3d Rz;

	Rz << cos(thetaz), -sin(thetaz), 0,
		sin(thetaz), cos(thetaz), 0,
		0, 0, 1;
	Matrix3d Rx;
	Rx <<	1, 0, 0,
			0, cos(thetax), -sin(thetax),
			0, sin(thetax), cos(thetax);
	R = Rx*Rz;
	this->p << 1.0, 1.0, 1.0;
	computeE();
	
	//w, v
	this->Omega.setZero();
	this->Omega << 0.0, 0.0, 0.0;
	this->V.setZero();
	this->V << 0.0, 0.0, 0.0;
	this->Phi.segment<3>(0) = Omega;
	this->Phi.segment<3>(3) = V;

	computePhiT();

}

void RigidBody::step(double h) {
	
	computeB(); // ->R->E

	shared_ptr<QuadProgMosek> program_ = make_shared <QuadProgMosek>();
	program = program_;
	program_->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
	program_->setParamInt(MSK_IPAR_LOG, 10);
	program_->setParamInt(MSK_IPAR_LOG_FILE, 1);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-8);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-10);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-8);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-8);
	program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-8);
	program->setNumberOfVariables(6);
	program->setLowerVariableBound(xl);
	program->setUpperVariableBound(xu);

	program->setObjectiveMatrix(A);

	RHS = M * Phi + h * (PhiT * M * Phi + B);
	program->setObjectiveVector(-RHS);
	bool success = program->solve();
	
	sol = program->getPrimalSolution();
	
	Phi = sol;
	Omega = Phi.segment<3>(0);
	V = Phi.segment<3>(3);
	computePhiT();
	computeEE();
	Etemp = E * (h * EE).exp();
	
	numCol = 0;

	// Update positions of all the nodes to detect collisions
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->xo = nodes[i]->x;
		nodes[i]->x = local2world(Etemp, nodes[i]->x0);
		if (nodes[i]->x(1) < -5.0) {
			index.push_back(i);
			cout << "col point:" << i << endl;
			numCol += 1;
		}
	}

	if (numCol != 0) {
		
		cout << "Phi: " << endl << Phi << endl;
		// Add constraints on the collision points
		program->setNumberOfInequalities(numCol);
		C.resize(numCol, 6);
		b.resize(numCol);
		b.setZero();

		VectorXd convec;
		convec.resize(6);
		//R = Etemp.block(0, 0, 3, 3);
		for (int c = 0; c < numCol; c++) {
			temp.block(0, 0, 3, 3) = vec2crossmatrix(nodes[index[c]]->x0).transpose();
			
			convec = ynormal.transpose() * R * temp;
			for (int id = 0; id < 6; id++) {
				C_.push_back(ETriplet(c, id, -convec(id)));
			}
		}
		C.setFromTriplets(C_.begin(), C_.end());
		program->setInequalityMatrix(C);
		program->setInequalityVector(b);

		C_.clear();
		index.clear();

		bool success = program->solve();
		sol = program->getPrimalSolution();
		Phi = sol;
		cout << "Phi:after " << endl << Phi << endl;
		
		Omega = Phi.segment<3>(0);
		V = Phi.segment<3>(3);
		computePhiT();
		computeEE();
		E = E * (h * EE).exp();
	}
	else {
		E = Etemp;
	}

	// E
	updateLocalFrame(); // update R, p
	updatePosNor(E);
}

void RigidBody::init() {
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &norBufID);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);

	glGenBuffers(1, &eleBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, eleBuf.size() * sizeof(unsigned int), &eleBuf[0], GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	assert(glGetError() == GL_NO_ERROR);
}

void RigidBody::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p)const {

	glUniform3fv(p->getUniform("kdFront"), 1, Vector3f(1.0, 0.0, 0.0).data());
	glUniform3fv(p->getUniform("kdBack"), 1, Vector3f(1.0, 1.0, 0.0).data());
	MV->pushMatrix();
	glUniformMatrix4fv(p->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	int h_pos = p->getAttribute("aPos");
	glEnableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);

	glBufferData(GL_ARRAY_BUFFER, posBuf.size() * sizeof(float), &posBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);

	int h_nor = p->getAttribute("aNor");
	glEnableVertexAttribArray(h_nor);
	glBindBuffer(GL_ARRAY_BUFFER, norBufID);
	glBufferData(GL_ARRAY_BUFFER, norBuf.size() * sizeof(float), &norBuf[0], GL_DYNAMIC_DRAW);
	glVertexAttribPointer(h_nor, 3, GL_FLOAT, GL_FALSE, 0, (const void *)0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, eleBufID);
	glDrawElements(GL_TRIANGLES, 3 * nTriFaces, GL_UNSIGNED_INT, (const void *)(0 * sizeof(unsigned int)));

	glDisableVertexAttribArray(h_nor);
	glDisableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	MV->popMatrix();
}

void RigidBody::updatePosNor(MatrixXd E) {
	// update 
	for (int i = 0; i < nodes.size(); i++) {
		nodes[i]->xo = nodes[i]->x;
		nodes[i]->x = local2world(E, nodes[i]->x0);
	}

	for (int iface = 0; iface < nTriFaces; iface++) {
		Vector3d p1 = nodes[out.trifacelist[3 * iface + 0]]->x;
		Vector3d p2 = nodes[out.trifacelist[3 * iface + 1]]->x;
		Vector3d p3 = nodes[out.trifacelist[3 * iface + 2]]->x;

		//Position
		Vector3d e1 = p2 - p1;
		Vector3d e2 = p3 - p1;
		Vector3d normal = e1.cross(e2);
		normal.normalize();

		for (int idx = 0; idx < 3; idx++) {
			posBuf[9 * iface + 0 + idx] = p1(idx);
			posBuf[9 * iface + 3 + idx] = p2(idx);
			posBuf[9 * iface + 6 + idx] = p3(idx);
			norBuf[9 * iface + 0 + idx] = normal(idx);
			norBuf[9 * iface + 3 + idx] = normal(idx);
			norBuf[9 * iface + 6 + idx] = normal(idx);
		}
	}
}

 Vector3d RigidBody::local2world(MatrixXd E, Eigen::Vector3d x) {

	 // homo coor
	 Eigen::VectorXd xh;
	 xh.resize(4);
	 xh.segment<3>(0) = x;
	 xh(3) = 1.0;

	 Eigen::Vector3d xw = (E * xh).segment<3>(0);
	 return xw;
}

 Matrix3d RigidBody::vec2crossmatrix(Vector3d a) {
	 Matrix3d A;
	 A << 0, -a(2), a(1),
		 a(2), 0, -a(0),
		 -a(1), a(0), 0;
	 return A;
 }

 void RigidBody::computeE() {
	 E.setZero();
	 E.block(0, 0, 3, 3) = R;
	 E.block(0, 3, 3, 1) = p;
	 E(3, 3) = 1.0;

 }

 void RigidBody::computeB() {
	 B.setZero();

	 // if only gravity is involved
	 Vector3d g;
	 g << 0, -9.8, 0;
	 B.segment<3>(3) = R.transpose() * mass * g;
 }

 void RigidBody::computeOMEGA() {
	OMEGA = vec2crossmatrix(Omega);
	
 }

 void RigidBody::computePHI() {
	 computeOMEGA();
	 PHI.block(0, 0, 3, 3) = OMEGA;
	 PHI.block(3, 3, 3, 3) = OMEGA;

 }

 void RigidBody::computePhiT(){ 
	 computePHI();
	 PhiT = PHI.transpose();

 }

 void RigidBody::computeEE() {
	 EE.block(0, 0, 3, 3) = OMEGA;
	 EE.block(0, 3, 3, 1) = V;
 }

 void RigidBody::updateLocalFrame() {
	 R = E.block(0, 0, 3, 3);
	 p = E.block(0, 3, 3, 1);
 }