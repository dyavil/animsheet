/*
 * CalculsMSS.cpp :
 * Copyright (C) 2016 Florence Zara, LIRIS
 *               florence.zara@liris.univ-lyon1.fr
 *               http://liris.cnrs.fr/florence.zara/
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/** \file CalculsMSS.cpp
Programme calculant pour chaque particule i d un MSS son etat au pas de temps suivant 
 (methode d 'Euler semi-implicite) : principales fonctions de calculs.
\brief Fonctions de calculs de la methode semi-implicite sur un systeme masses-ressorts.
*/ 

#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "vec.h"
#include "ObjetSimule.h"
#include "ObjetSimuleMSS.h"
#include "Viewer.h"

using namespace std;





/**
* Calcul des forces appliquees sur les particules du systeme masses-ressorts.
 */
void ObjetSimuleMSS::CalculForceSpring()
{
	/// f = somme_i (ki * (l(i,j)-l_0(i,j)) * uij ) + (nuij * (vi - vj) * uij) + (m*g) + force_ext
	
	/// Rq : Les forces dues a la gravite et au vent sont ajoutees lors du calcul de l acceleration
    for (int i = 0; i < _SytemeMasseRessort->GetNbParticule(); ++i)
    {
        Particule *currentPart = _SytemeMasseRessort->GetPartList()[i];
    	Vector somme_i = Vector(0.0,0.0,0.0);
    	for(int j=0; j < currentPart->GetNbVoisins(); ++j){
            
            Ressort *r = currentPart->GetRessortList()[j];
            if(currentPart->GetRessortList().size() == 1 ) std::cout << "taille : " << currentPart->GetRessortList().size() <<std::endl;

            Particule *partA = r->GetParticuleA();
            Particule *partB = r->GetParticuleB();

            //si la particule courante est B on inverse
            if(r->GetParticuleB()->_Id == currentPart->_Id){
                partA = r->GetParticuleB();
                partB = r->GetParticuleA();
            }

    

            //lg ressort
    		float distj = distance(Point(P[partA->_Id]), Point(P[partB->_Id]));
            if(distj > 0.85) {
                //std::cout << "dechire : " << distj << endl;
                //currentPart->GetRessortList().erase(currentPart->GetRessortList().begin()+j);
                //j--;
            }

            //ki * (l(i,j)-l_0(i,j)))
    		float scal = (r->GetSpring()->_Raideur*(distj - r->GetSpring()->_L0));

            //uij
    		Vector elast = scal * normalize(P[partB->_Id]-P[partA->_Id]);//(P[partB->_Id]-P[partA->_Id])/distj;


    		float amorti = dot(Vector(r->GetSpring()->_Amorti, r->GetSpring()->_Amorti, r->GetSpring()->_Amorti), 
    			Vector((P[partB->_Id].x-P[partA->_Id].x)/distj, 
    			(P[partB->_Id].y-P[partA->_Id].y)/distj, 
    			(P[partB->_Id].z-P[partA->_Id].z)/distj)); 

    		Vector visc = Vector((P[partB->_Id].x-P[partA->_Id].x)/distj*amorti, 
    			(P[partB->_Id].y-P[partA->_Id].y)/distj*amorti, 
    			(P[partB->_Id].z-P[partA->_Id].z)/distj*amorti);
    		//std::cout << "ddddd " << elast.x << "   " << scal << std::endl;
    		elast = elast + visc;
    		
            /*if (elast.x+elast.y+elast.z > 12000)
            {
                std::cout << "ggggggg " << elast.x+elast.y+elast.z << "   " << std::endl;
                currentPart->GetRessortList().erase(currentPart->GetRessortList().begin()+j-1);
            }*/
            somme_i = somme_i + elast;
            
            //std::cout << "in " << somme_i << "   " << distj << "   " << scal << std::endl;
    	}
        Force[currentPart->_Id] = somme_i;
    	//Force.assign(i, somme_i);
        //std::cout << "ggggggg " << somme_i << "   " << std::endl;
    	
    }
		
}//void


/**
 * Gestion des collisions avec le sol - plan (x,y,z).
 */
void ObjetSimuleMSS::CollisionPlan(float x, float y, float z)
{
    /// Arret de la vitesse quand touche le plan
   for (int i = 0; i < _SytemeMasseRessort->GetNbParticule(); ++i){
        Particule *currentPart = _SytemeMasseRessort->GetPartList()[i];
        if (P[currentPart->_Id].x <= x || P[currentPart->_Id].y <= y+0.1 || P[currentPart->_Id].z <= z)
        {
            if(V[currentPart->_Id].y != 0) V[currentPart->_Id] = Vector(0, 0, 0);
            //std::cout << V[currentPart->_Id] << std::endl;
        }
   }
   
    
}// void

void ObjetSimuleMSS::CollisionTable(Point pmin, Point pmax){
    for (int i = 0; i < _SytemeMasseRessort->GetNbParticule(); ++i){
        Particule *currentPart = _SytemeMasseRessort->GetPartList()[i];
        if (((P[currentPart->_Id].x >= pmin.x && P[currentPart->_Id].x <= pmax.x) && ( P[currentPart->_Id].z >= pmin.z && P[currentPart->_Id].z <= pmax.z)) && P[currentPart->_Id].y <= pmin.y)
        {
            //if (P[currentPart->_Id].z <= -1) std::cout << P[currentPart->_Id] << std::endl;
            if(V[currentPart->_Id].y != 0) V[currentPart->_Id] = Vector(0, 0, 0);
        }
        //std::cout << "down" << std::endl;
   }
}
