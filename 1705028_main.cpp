#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<string>
#include <windows.h>
#include <GL/glut.h>
#include "header.h"
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))
#define theta (3.0/180.0)
using namespace std;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

double dis;
int slices;
int stacks;
int floorWidth =1000;
int tileWidth =20;
int recursion;
int pixels;
int objects_count;
int window_width = 500;
int window_height = 500;
double view_angle = 80.0;
int image_count=0;
extern Point pos;
extern int recur_level;
extern vector <Object*> objects;
extern vector <PointLight> pointLights;
extern vector <SpotLight> spotLights;


struct point
{
	double x,y,z;
};

struct point u;
struct point r;
struct point l;

class cross_vec
{
    struct point product;
    struct point a;
    struct point b;
public:
    cross_vec(struct point p, struct point q )
    {
        a = p;
        b = q;
    }
    struct point operate ()
    {
        product.x = a.y * b.z - a.z * b.y;
        product.y = a.z * b.x - a.x * b.z;
        product.z = a.x * b.y - a.y * b.x;
        return product;
    }

};
void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 200,0,0);
			glVertex3f(-200,0,0);

			glVertex3f(0,-200,0);
			glVertex3f(0, 200,0);

			glVertex3f(0,0, 200);
			glVertex3f(0,0,-200);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}


void capture() {
    cout << pos.x<<" "<<pos.y<<" "<<pos.z<< " : Start" << endl;
    bitmap_image bit_img(pixels, pixels);

    for(int i=0; i<pixels; i++) ///column
    {
        for(int j=0; j<pixels; j++) ///row
         {
            bit_img.set_pixel(i, j, 0.0, 0.0, 0.0);  // background color = black
        }
    }
    double planeDistance = (window_height/(2.0))/(tan(view_angle/2.0*(pi/180.0)));
    double width = (window_width/2.0);
    double height = (window_height/2.0);
    Point tempL(l.x,l.y,l.z);
    Point tempR(r.x,r.y,r.z);
    Point tempU(u.x,u.y,u.z);
    Point topLeft = pos+ tempL*planeDistance-tempR*width+tempU*height;


    double du = ((double) window_width/pixels);
    double dv = ((double) window_height/pixels);

    /* middle of the grid cell */
    topLeft = topLeft+ tempR*(du/2.0)-tempU*(dv/2.0);


    for(int col=0; col<pixels; col++) {
        for(int row=0; row<pixels; row++) {
             //calculating current pixel and casting ray from camera to (curPixel-camera) direction
            Point curPixel = topLeft+tempR*(col*du)-tempU*(row*dv);
            Ray ray(pos, curPixel-pos); // camera to (curPixel-camera) direction

             //finding nearest intersecting object (if available)
            int near_idx = INT_MAX;
            double t, tMin=INF;

            for(int i=0; i<objects.size(); i++) {
                Color color;  // dummyColor
                t = objects[i]->intersect(ray,0,color);

                if(t>0.0 && t<tMin) {
                    tMin = t;
                    near_idx = i;
                }
            }
            if(near_idx != INT_MAX) {
                //pixel print ekhane kore
                Color color;
                tMin = objects[near_idx]->intersect(ray, 1, color);
                int r,g,b;
                r = (int) round(color.getR()*255.0);
                g = (int) round(color.getG()*255.0);
                b=(int) round(color.getB()*255.0);
                bit_img.set_pixel(col, row, r, g, b);
            }
        }
    }
    image_count+=1;
    bit_img.save_image("E:\\4-1\\OUR DRIVE\\Sessional\\CSE 310 computer graphics\\offline 3\\code_index_2\\output"+to_string(image_count)+".bmp");
    cout << pos.getX()<<" "<<pos.getY()<<" "<<pos.getZ()<<" "<< ": End" << endl;

}
void keyboardListener(unsigned char key, int x,int y){
	switch(key){
        case '0':
            capture();
		case '1':
            {
            cross_vec temp(u,r);
            cross_vec temp2(l,u);
		    r.x = r.x*cos(theta) +  temp.operate().x* sin(theta) ;
		    r.y = r.y*cos(theta) +  temp.operate().y* sin(theta) ;
		    r.z = r.z*cos(theta) +  temp.operate().z* sin(theta) ;

		    l.x = l.x* cos(theta) - temp2.operate().x*sin(theta);
            l.y = l.y* cos(theta) - temp2.operate().y*sin(theta) ;
            l.z = l.z* cos(theta) - temp2.operate().z*sin(theta) ;
			break;
            }
        case '2':
			{
            cross_vec temp(u,r);
            cross_vec temp2(l,u);
		    r.x = r.x*cos(theta) -  temp.operate().x* sin(theta) ;
		    r.y = r.y*cos(theta) -  temp.operate().y* sin(theta) ;
		    r.z = r.z*cos(theta) -  temp.operate().z* sin(theta) ;

		    l.x = l.x* cos(theta) + temp2.operate().x*sin(theta);
            l.y = l.y* cos(theta) + temp2.operate().y*sin(theta) ;
            l.z = l.z* cos(theta) + temp2.operate().z*sin(theta) ;
			break;
            }
        case '3':
            {
            cross_vec temp(u,r);
            cross_vec temp2(r,l);
		    u.x = u.x*cos(theta) -  temp.operate().x* sin(theta) ;
		    u.y = u.y*cos(theta) -  temp.operate().y* sin(theta) ;
		    u.z = u.z*cos(theta) -  temp.operate().z* sin(theta) ;

		    l.x = l.x* cos(theta) + temp2.operate().x*sin(theta);
            l.y = l.y* cos(theta) + temp2.operate().y*sin(theta) ;
            l.z = l.z* cos(theta) + temp2.operate().z*sin(theta) ;
            break;
            }
        case '4':
			{
            cross_vec temp(u,r);
            cross_vec temp2(r,l);
		    u.x = u.x*cos(theta) +  temp.operate().x* sin(theta) ;
		    u.y = u.y*cos(theta) +  temp.operate().y* sin(theta) ;
		    u.z = u.z*cos(theta) +  temp.operate().z* sin(theta) ;

		    l.x = l.x* cos(theta) - temp2.operate().x*sin(theta);
            l.y = l.y* cos(theta) - temp2.operate().y*sin(theta) ;
            l.z = l.z* cos(theta) - temp2.operate().z*sin(theta) ;
            break;
            }
        case '5':
			{
            cross_vec temp(r,l);
            cross_vec temp2(l,u);
		    r.x = r.x*cos(theta) -  temp.operate().x* sin(theta) ;
		    r.y = r.y*cos(theta) -  temp.operate().y* sin(theta) ;
		    r.z = r.z*cos(theta) -  temp.operate().z* sin(theta) ;

		    u.x = u.x* cos(theta) + temp2.operate().x*sin(theta);
            u.y = u.y* cos(theta) + temp2.operate().y*sin(theta) ;
            u.z = u.z* cos(theta) + temp2.operate().z*sin(theta) ;
            break;
            }
        case '6':
            {
            cross_vec temp(r,l);
            cross_vec temp2(l,u);
		    r.x = r.x*cos(theta) +  temp.operate().x* sin(theta) ;
		    r.y = r.y*cos(theta) +  temp.operate().y* sin(theta) ;
		    r.z = r.z*cos(theta) +  temp.operate().z* sin(theta) ;

		    u.x = u.x* cos(theta) - temp2.operate().x*sin(theta);
            u.y = u.y* cos(theta) - temp2.operate().y*sin(theta) ;
            u.z = u.z* cos(theta) - temp2.operate().z*sin(theta) ;
            break;
            }
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			pos.x = pos.x - l.x;
			pos.y = pos.y - l.y;
			pos.z = pos.z - l.z;
			break;
		case GLUT_KEY_UP:		// up arrow key

            pos.x = pos.x + l.x;
			pos.y = pos.y + l.y;
			pos.z = pos.z + l.z;
			break;
		case GLUT_KEY_RIGHT:
			pos.x = pos.x + r.x;
			pos.y = pos.y + r.y;
			pos.z = pos.z + r.z;
			break;
		case GLUT_KEY_LEFT:
			pos.x = pos.x - r.x;
			pos.y = pos.y - r.y;
			pos.z = pos.z - r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x = pos.x + u.x;
			pos.y = pos.y + u.y;
			pos.z = pos.z + u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos.x = pos.x - u.x;
			pos.y = pos.y - u.y;
			pos.z = pos.z - u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
		     {
//               increment_radius();
		       break;
		     }
			break;
		case GLUT_KEY_END:
		     {
//               decrement_radius();
		       break;
		     }


		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
//				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}


void loadData()
{
    ifstream input;
    input.open("E:\\4-1\\OUR DRIVE\\Sessional\\CSE 310 computer graphics\\offline 3\\code_index_2\\scene.txt");
    // populating in the loadData() function
    if(!input.is_open()) {
        cout << "input.is_open(): failed to open input file" << endl;
        exit(EXIT_FAILURE);
    }
    input>>recur_level;
    input>>pixels;
    input>>objects_count;
    int check=0;
    Object *obj=NULL;
    cout<<objects_count<<endl;
    for(int i=0; i<objects_count; i++)
    {
        string shape="";
        input>>shape;
        double x,y,z,r,g,b, ambient, diffuse, specular, reflection_coefficient, shininess;
         Color clr;
        if(shape == "sphere"){
            double radius;
            input>>x>>y>>z;
//            cout<<x<<" "<<y<<" "<<z<<endl;
            Point center(x,y,z);
            input>>radius;
            input>>x>>y>>z;
            input>>ambient>>diffuse>>specular>>reflection_coefficient;
            input>>shininess;

            clr.setColor(x,y,z);
            obj = new Sphere(center, radius,slices,stacks);
            obj->setColor(clr);
//            obj->setSliceAndStack();
            obj->setCoEfficients(ambient,diffuse,specular,reflection_coefficient);
            obj->setShine(shininess);
            objects.push_back(obj);
        }
        else if(shape == "triangle"){
            input>>x>>y>>z;
            Point a(x,y,z);
            input>>x>>y>>z;
            Point b(x,y,z);

            input>>x>>y>>z;
            Point c(x,y,z);

            double red,green,blue;
            input>>red>>green>>blue;
            input>>ambient>>diffuse>>specular>>reflection_coefficient;
            input>>shininess;

            clr.setColor(red,green,blue);
            obj = new Triangle(a,b,c);
            obj->setColor(clr);
            obj->setCoEfficients(ambient,diffuse,specular,reflection_coefficient);
            obj->setShine(shininess);
            objects.push_back(obj);
        }
        else if(shape == "general"){
            double A, B, C, D, E, F, G, H, I, J;
            double cube_ref_point, length, width, height;
            input>>A>>B>>C>>D>>E>>F>>G>>H>>I>>J;
            struct GeneralPoints p;
            p.a=A;p.b=B;p.c=C;p.d=D;p.e=E;p.f=F;p.g=G;p.h=H;p.i=I;p.j=J;
            input>>x>>y>>z;
            Point ref_p(x,y,z);
            input>>length>>width>>height;
            input>>r>>g>>b;
            input>>ambient>>diffuse>>specular>>reflection_coefficient;
            input>>shininess;

            clr.setColor(r,g,b);
            obj = new General(p,ref_p, length,width,height);
            obj->setColor(clr);
            obj->setCoEfficients(ambient,diffuse,specular,reflection_coefficient);
            obj->setShine(shininess);
            objects.push_back(obj);
        }
        else{
            cout <<shape<<" : object is invalid" << endl;
            check = 1;
            break;
        }

    }
//    if(check ==1) exit();
    ///Floor
    Color c;
    obj = new Floor(1000.0, 20.0,c);
    obj->setColor(Color(1.0, 1.0, 1.0));
    obj->setCoEfficients(0.25, 0.25, 0.25, 0.25);
    obj->setShine(15);
    objects.push_back(obj);
    obj=NULL;

    int point_light_count=0;
    input>>point_light_count;
    PointLight pl;
    for(int i=0; i< point_light_count; i++){
            double x,y,z, red,green,blue;
            Color clr;Point p;

            input>>x>>y>>z; // taking position
            input>>red>>green>>blue;//taking color

            p.setPoint(x,y,z);
            clr.setColor(red,green,blue);
            double radius=1.0;

            pointLights.push_back(PointLight(p,clr,radius, slices, stacks)); //push into object
    }
    int spot_light_count=0;
    input>>spot_light_count;
    cout<<spot_light_count<<endl;
    for(int i=0; i< spot_light_count; i++){
            PointLight pl;
            double x,y,z, red,green,blue, dirx,diry,dirz,cutoff_angle;
            Color clr;
            Point p, p_dir;

            input>>x>>y>>z;
            input>>red>>green>>blue;
            input>>dirx>>diry>>dirz;
            input>>cutoff_angle;

            p.setPoint(x,y,z);
            clr.setColor(red,green,blue);
            p_dir.setPoint(dirx,diry,dirz);
            pl.setPointAndColor(p,clr);
            double radius = 10;
            spotLights.push_back(SpotLight(pl,p_dir, cutoff_angle,radius,slices,stacks));
    }
    input.close();


}
void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);
	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
//	gluLookAt(50*cos(cameraAngle), 50*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);
	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);
	/****************************
	/ Add your objects from here
	****************************/
	drawGrid();

    for(int i=0; i<objects.size(); i++) {
        objects[i]->draw();
	}

	for(int i=0; i<pointLights.size(); i++) {
        pointLights[i].draw();
	}
	for(int i=0; i<spotLights.size(); i++) {
        spotLights[i].draw();
	}
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}


void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=0.0;
	cameraAngle=1.0;
	angle=0;


    slices = 24;
    stacks = 30;
    Object *obj;

    objects_count=0;

	//clear the screen
	glClearColor(0,0,0,0);

	///************************
	pos.x = 100.0;
	pos.y = 100.0;
	pos.z = 30.0;

	u.x = 0;
	u.y = 0;
	u.z = 1.0;

	r.x = -1/sqrt(2);
	r.y = 1/sqrt(2);
	r.z = 0;

	l.x = -1/sqrt(2);
	l.y = -1/sqrt(2);
	l.z = 0;
	///************************
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
}
void lightMemoryFree()
{
    for(int i=0; i<objects.size(); i++) {
        delete objects[i];
    }
    objects.clear();
}
void objectMemoryFree()
{
    pointLights.clear();
    spotLights.clear();
}
int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	loadData();
	int value1, value2;
    value1 = atexit(lightMemoryFree);
    value2 = atexit(objectMemoryFree);

	if((value1 != 0) || (value2 != 0)) {
        cout << "Programme falied" << endl;
        exit(EXIT_FAILURE);
	}

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
