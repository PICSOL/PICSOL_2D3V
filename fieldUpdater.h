#pragma once
#include "Global.h"
#include "scalarField.h"
#include "vectorField.h"
#include "Vector3D.h"
#include "Species.h"
#include "ellipticSolver.h"
#include "fieldOperator.h"

/*
  This field updater is the most essential part of PICSOL code, which updates the electromagnetic field using 
  Vlasov-Darwin equations, all physical quantities have been normalized, and boundary condition for the eletromagnetic 
  field is dirichlet in x direction (index is i), periodic in y direction (index is j).
   
    field that requires update : rho, J, phi ,E_l, E_t, E_total, B

  Edited by Huaxiang Zhang 2019/4/3
*/

/*
rhoComp - compute total rho by adding up electron and ion contribution
*/
void rhoComp(scalarField& electDen, scalarField& ionDen, scalarField& rho);

/*
JComp - compute total current density J by adding up electron and ion contribution
*/
void JComp(vectorField& electCur, vectorField& ionCur, vectorField& current);

/*
JComp - compute total current density J by adding up electron and ion contribution
*/
void convectComp(vectorField& electCon, vectorField& ionCon, vectorField& convect);

/*
EtotalComp - compute total electric field by adding up longtitude part and transverse part
*/
void EtotalComp(vectorField& E_l, vectorField& E_t, vectorField& E_total);

/*
phiComp - compute phi field by solving one 2D poisson equation
*/
void phiComp(scalarField &rho, scalarField &phi);

/*
ElComp - compute E_l field by computing the divergence of phi
*/
void ElComp(scalarField &phi, vectorField &E_l);

/*
BComp - compute B field by solving three 2D-poisson equations, B = B_0 + B_1,
        we only compute B_1, B_0 is the background static magnetic field 
		much greater than B_1
*/
void BComp(vectorField &J, vectorField &B);

/*
EtComp - compute E_t field by solving three 2D-helmholtz equations, and one poisson equation
*/
void EtComp(scalarField &rho, vectorField &Current, vectorField &Convect, 
	      vectorField &E_l, vectorField &B, vectorField &E_t);



