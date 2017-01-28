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
    	Vector somme_i = Vector(0.0,0.0,0.0);
    	for(int j=0; j < _SytemeMasseRessort->GetPartList()[i]->GetNbVoisins(); ++j){
    		float distj = sqrt(pow((_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().x-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().x), 2)
    			+pow((_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().y-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().y), 2)
    			+pow((_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().z-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().z), 2));
    		float scal = (_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetSpring()->_Raideur*(distj
    			-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetSpring()->_L0));
    		Vector elast = Vector((_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().x-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().x)/distj*scal, 
    			(_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().y-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().y)/distj*scal, 
    			(_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().z-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().z)/distj*scal);

    		float amorti = dot(Vector(_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetSpring()->_Amorti, _SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetSpring()->_Amorti, _SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetSpring()->_Amorti), 
    			Vector((_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().x-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().x)/distj, 
    			(_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().y-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().y)/distj, 
    			(_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().z-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().z)/distj)); 

    		Vector visc = Vector((_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().x-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().x)/distj*amorti, 
    			(_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().y-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().y)/distj*amorti, 
    			(_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleB()->GetPosition().z-_SytemeMasseRessort->GetPartList()[i]->GetRessortList()[j]->GetParticuleA()->GetPosition().z)/distj*amorti);
    		//std::cout << "ddddd " << elast.x << "   " << scal << std::endl;
    		elast = elast + visc;
    		somme_i = somme_i + elast;

    	}
    	Force.assign(i, somme_i);
    	//std::cout << "ggggggg " << somme_i.x << "   " << _VISize << std::endl;
    }
		
}//void


/**
 * Gestion des collisions avec le sol - plan (x,y,z).
 */
void ObjetSimuleMSS::CollisionPlan(float x, float y, float z)
{
    /// Arret de la vitesse quand touche le plan
   
    
}// void

