/** \file Viewer-init.cpp
 * \brief Fonctions d initialisation des Mesh des objets a afficher.
 */

#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "draw.h"
#include "wavefront.h"
#include "texture.h"
#include "Viewer.h"

using namespace std;





/*
 * Creation du maillage des axes pour les visualiser facilement.
 */
void Viewer::init_axe()
{
    m_axe = Mesh(GL_LINES);
    
    // Axe x en rouge
    m_axe.color( Color(1, 0, 0));
    m_axe.vertex( 0,  0, 0);
    m_axe.vertex( 1,  0, 0);
    
    // Axe y en vert
    m_axe.color( Color(0, 1, 0));
    m_axe.vertex( 0,  0, 0);
    m_axe.vertex( 0,  1, 0);
    
    // Axe z en bleu
    m_axe.color( Color( 0, 0, 1));
    m_axe.vertex( 0,  0, 0);
    m_axe.vertex( 0,  0, 1);
    
}


/*
 * Creation du maillage d'une grille.
 */
void Viewer::init_grid()
{
    m_grid = Mesh(GL_LINES);
    
    // Couleur de la grille
    m_grid.color( Color(1, 1, 1));
    
    for(int i=-5;i<=5;++i)
        for(int j=-5;j<=5;++j)
        {
            m_grid.vertex( -5, 0, j);
            m_grid.vertex( 5, 0,  j);
            
            m_grid.vertex( i, 0, -5);
            m_grid.vertex( i, 0, 5);
            
        }
}


/*
 * Creation du maillage du plan de collision.
 */
void Viewer::init_plan(float x, float y, float z)
{
    m_plan = Mesh(GL_TRIANGLE_STRIP);
    
    // Couleur de la grille
    //m_plan.color( Color(1, 1, 0));
    
    /*for(int i=-5;i<=5;++i)
        for(int j=-5;j<=5;++j)
        {
            m_plan.vertex( -5+x, y, j+z);
            m_plan.vertex( 5+x, y,  j+z);
            
            m_plan.vertex( i+x, y, -5+z);
            m_plan.vertex( i+x, y, 5+z);
            
        }*/
    float t = 0.1;
    for (int j = -5; j < 5; ++j)
    {
        
        for (int i = -5; i < 5; ++i)
        {
            m_plan.texcoord(0, 0);
            int a = m_plan.vertex(i, y, j+1);
            m_plan.texcoord(1, 1);
            int b = m_plan.vertex(i+1, y, j);
            m_plan.texcoord(1, 0);
            int c = m_plan.vertex(i, y, j);
            m_plan.texcoord(0, 1);
            int d = m_plan.vertex(i+1, y, j+1);
            m_plan.triangle(a, c, b);
            m_plan.triangle(a, b, d);
        }
        if(t > 0) t = 0;
        else t = 0.1;
        m_plan.restart_strip();
    }
    /*m_plan.texcoord(0, 0);
    int a = m_plan.vertex(-5, y, 5);
    m_plan.texcoord(1, 1);
    int b = m_plan.vertex(5, y, -5);
    m_plan.texcoord(1, 0);
    int c = m_plan.vertex(-5, y, -5);
    m_plan.texcoord(0, 1);
    int d = m_plan.vertex(5, y, 5);
    m_plan.triangle(a, c, b);
    m_plan.triangle(a, b, d);*/
}

void Viewer::init_scene()
{
    Mesh m1 = read_mesh("/home/dyavil/Documents/Master/SyntheseImage/gkit2light/data/bigguy.obj");
    //m1.color(1, 0, 0);
    Mesh m2 = read_mesh("/home/dyavil/Documents/Master/SyntheseImage/gkit2light/data/bigguy.obj");
    Mesh m3 = read_mesh("data/Barrel.obj");
    m_scene.push_back(m1);
    m_scene.push_back(m2);
    m_scene.push_back(m3);
}

/*
 * Creation du maillage d un cube.
 */
void Viewer::init_cube()
{
    // Definition des 8 sommets du cube
    //                          0           1           2       3           4           5       6           7
    // 8 sommets - chaque sommet a coodonnees en x, y, z
    static float pt[8][3] = { {-1,-1,-1}, {1,-1,-1}, {1,-1,1}, {-1,-1,1}, {-1,1,-1}, {1,1,-1}, {1,1,1}, {-1,1,1} };
   
    // 6 faces - chacune des faces a 4 sommets
    static int f[6][4] = {    {0,1,2,3}, {5,4,7,6}, {2,1,5,6}, {0,3,7,4}, {3,2,6,7}, {1,0,4,5} };
    
    // 6 normales (une normale par face) - chaque normale a coodonnees en x, y, z
    static float n[6][3] = { {0,-1,0}, {0,1,0}, {1,0,0}, {-1,0,0}, {0,0,1}, {0,0,-1} };

    
    // Maillage de type GL_TRIANGLE_STRIP
    m_cube = Mesh(GL_TRIANGLE_STRIP);
    
    // Couleur du cube
    m_cube.color( Color(1, 1, 1) );
    
    // Boucle sur les 6 faces du cube
    for (int i=0;i<6;i++)
    {
        // Definition de la normale a la face i
        m_cube.normal(  n[i][0], n[i][1], n[i][2] );
        
        // Definition du sommet 0 de la face i : coordonnees texture + coordonnees geometrique
        m_cube.texcoord( 0,0 );
        m_cube.vertex( pt[ f[i][0] ][0], pt[ f[i][0] ][1], pt[ f[i][0] ][2] );
        
        // Definition du sommet 1 de la face i  : coordonnees texture + coordonnees geometrique
        m_cube.texcoord( 1,0);
        m_cube.vertex( pt[ f[i][1] ][0], pt[ f[i][1] ][1], pt[ f[i][1] ][2] );
        
        // Definition du sommet 3 de la face i  : coordonnees texture + coordonnees geometrique
        m_cube.texcoord(0,1);
        m_cube.vertex(pt[ f[i][3] ][0], pt[ f[i][3] ][1], pt[ f[i][3] ][2] );
        
        // Definition du sommet 2 de la face i : coordonnees texture + coordonnees geometrique
        m_cube.texcoord(1,1);
        m_cube.vertex( pt[ f[i][2] ][0], pt[ f[i][2] ][1], pt[ f[i][2] ][2] );
        
        m_cube.restart_strip();
        
    }
}

void Viewer::init_table()
{
    m_table = read_mesh("data/table.obj");
    /*m_table = Mesh(GL_TRIANGLE_STRIP);
    m_table.color( Color(0.4, 0.2, 0.8) );
    m_table.texcoord( 0,0 );
    int a = m_table.vertex(Point(0.39, 0.82, 0.5));
    m_table.texcoord( 1,1 );
    int b = m_table.vertex(Point(1.42, 0.82, -0.5));
    m_table.texcoord( 1,0 );
    int c = m_table.vertex(Point(0.39, 0.82, -0.5));
    m_table.texcoord( 0,1 );
    int d = m_table.vertex(Point(1.42, 0.82, 0.5));
    m_table.triangle(a, c, b);
    m_table.triangle(a, b, d);

    m_table.restart_strip();
    m_table.texcoord( 0,0 );
    int e = m_table.vertex(Point(0.39, 0.35, 0.5));
    m_table.texcoord( 0,1 );
    int f = m_table.vertex(Point(1.42, 0.35, -0.5));
    m_table.texcoord( 0,0 );
    int g = m_table.vertex(Point(0.39, 0.35, -0.5));
    m_table.texcoord( 0,1 );
    int h = m_table.vertex(Point(1.42, 0.35, 0.5));
    m_table.triangle(e, f, g);
    m_table.triangle(e, h, f);

    m_table.restart_strip();
    m_table.triangle(e, a, d);
    m_table.triangle(e, d, h);
    m_table.restart_strip();
    m_table.triangle(h, d, b);
    m_table.triangle(h, b, f);
    m_table.restart_strip();
    m_table.triangle(f, b, c);
    m_table.triangle(f, c, g);
    m_table.restart_strip();
    m_table.triangle(g, c, a);
    m_table.triangle(g, a, e);*/
    //Point t1, t2;
    /*m_table.bounds(t1, t2);
    std::cout << t1 << t2 << std::endl;*/
    m_table_pos = make_pair(Point(1, -7, -2), Point(5, -7, 2));
}


/*
 * Creation du maillage d'une sphere - centre (0, 0, 0) - rayon = r.
 */
void Viewer::init_sphere()
{
    const int divBeta = 26;
    const int divAlpha = divBeta/2;
    int i,j;
    float beta, alpha, alpha2;
    float rayon = 0.05;
    
    m_sphere = Mesh(GL_TRIANGLE_STRIP);
    
    
    m_sphere.color( Color(0, 0, 1) );
            
    for(i=0;i<divAlpha;++i)
    {
        alpha = -0.5f*M_PI + float(i)*M_PI/divAlpha;
        alpha2 = -0.5f*M_PI + float(i+1)*M_PI/divAlpha;
        
        for(j=0;j<divBeta;++j)
        {
            beta = float(j)*2.f*M_PI/(divBeta-1);
            
            m_sphere.texcoord(beta/(2.0f*M_PI), 0.5f+alpha2/M_PI);
            
            m_sphere.normal( - Vector(rayon * cos(alpha2)*cos(beta),
                                    rayon * sin(alpha2),
                                    rayon * cos(alpha2)*sin(beta)) );
            
            m_sphere.vertex( Point(rayon * cos(alpha2)*cos(beta),
                                   rayon*sin(alpha2),
                                   rayon * cos(alpha2)*sin(beta))   );
            
            
            m_sphere.texcoord(beta/(2.0f*M_PI), 0.5f+alpha/M_PI);
            
            m_sphere.normal( - Vector(rayon * cos(alpha)*cos(beta),
                                    rayon * sin(alpha),
                                    rayon * cos(alpha)*sin(beta)) );
            
            m_sphere.vertex( Point(rayon * cos(alpha)*cos(beta),
                                   rayon * sin(alpha),
                                   rayon * cos(alpha)*sin(beta)) );
        }
        
        m_sphere.restart_strip();
    }
    
    
}


/*
 * Fonction d initialisation.
 */
int Viewer::init()
{
    // Etats par defaut openGL
    glClearColor(0.5f, 0.6f, 0.8f, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Camera
    m_camera.lookat( Point(0,0,0), 15 );
    // m_camera.lookat( Point(0,0,0), 150 );
    
    // Lumiere
    gl.light( Point(0, 20, 20), White() );
    
    // Chargement des textures utilisees dans la scene
    // Texture pour le cube
   //  m_cube_texture = read_texture(0, "../data/debug2x2red.png");
    m_cube_texture = read_texture(0, smart_path("data/debug2x2red.png"));

    // Texture pour le tissu
    //m_tissu_texture = read_texture(0, "data/papillon.png");
    //m_tissu_texture = read_texture(0, "data/textures/tissu1.png");
   // m_tissu_texture = read_texture(0, "data/textures/tissu2.jpg");
    m_tissu_texture = read_texture(0, smart_path("data/textures/tissu1.png"));
    m_plan_texture = read_texture(0, "data/floor1.jpg");
    m_bg_texture = read_texture(0, "data/fire.png");
    m_table_texture = read_texture(0, "data/floor1.jpg");
    m_obj_texture1 = read_texture(0, "data/wood.jpg");
    m_obj_texture2 = read_texture(0, "data/liege.jpg");
    m_obj_texture3 = read_texture(0, "data/metal.jpg");

    // Appel des procedures d initialisation des objets de la scene
    // Pour les objets non simules
    // Exemple : init_votreObjet();
    init_axe();
    init_grid();
    init_cube();
    init_sphere();
    init_table();
    init_scene();
    
    // Creation du plan (x, y, z) - plan utilise pour les ObjetSimule::CollisionPlan(x, y, z);
    // Rq : pas vraiment le plan, mais < x, < y, < z
    init_plan(0, 0, 0);
    
    // Initialisation du Tps
    Tps = 0;
    
    // Point interaction
    MousePos = Vector(0, 0, 0);
    
    return 0;
    
}


