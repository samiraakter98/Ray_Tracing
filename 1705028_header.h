#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<limits>

#include <windows.h>
#include <GL/glut.h>

using namespace std;
#define pi (2*acos(0.0))
#define INF numeric_limits<double>::infinity()

int recur_level=0;
struct GeneralPoints
{
    double a,b,c,d,e,f,g,h,i,j;
};
class Color{
    double r,g,b;
public:
    Color(){r=g=b=0.0;}
    Color(double r,double g, double b){
        this->r=r;
        this->g=g;
        this->b=b;
    }
    void setColor(double r,double g, double b){
        this->r=r;
        this->g=g;
        this->b=b;
    }
    double getR(){return r;}
    double getG(){return g;}
    double getB(){return b;}
};

class Point{
public:
    double x,y,z;
    Point(){
        x=0.0;
        y=0.0;
        z=0.0;
    }
    Point(double x, double y, double z){
        this->x = x;this->y=y;this->z=z;
    }
    double getX(){return x;}
    double getY(){return y;}
    double getZ(){return z;}
    void setX(double x){ this->x =x;}
    void setY(double y){this->y=y;}
    void setZ(double z){this->z=z;}
    void setPoint(double x, double y, double z){ this->x = x;this->y=y;this->z=z;}
    double dotMultiplication(Point p) {
        return x*p.getX() + y*p.getY() + z*p.getZ();
    }
    Point crossMultiplication(Point p) {
        Point temp;
        temp.setX(y * p.getZ() - z * p.getY());
        temp.setY(z* p.getX() - x * p.getZ()) ;
        temp.setZ(x* p.getY() - y * p.getX());
        return temp;
    }
    Point operator+(Point p) {
        Point temp;
        temp.setPoint(x + p.getX() ,  y + p.getY(), z + p.getZ());
        return temp;
    }
    Point operator-(Point p) {
        Point temp;
        temp.setPoint(x - p.getX() ,  y - p.getY(), z - p.getZ());
        return temp;
    }
    Point operator*(double num) {
        Point temp;
        temp.setPoint(x*num , y*num, z*num);
        return temp;
    }
    double distance(Point p){
        double a = pow(x-p.getX(), 2);
        double b = pow(y-p.getY(), 2);
        double c = pow(z-p.getZ(), 2);
        return sqrt(a+b+c);
    }
    Point normalize(){
        double num = sqrt(pow(x, 2.0)+pow(y, 2.0)+pow(z, 2.0));;
        x = x/num;
        y = y/num;
        z = z/num;
    }
    double getAngle(Point p)
    {
        double dot_value = x*p.getX() + y*p.getY() + z*p.getZ();

        double a = sqrt(pow(x, 2.0)+pow(y, 2.0)+pow(z, 2.0));
        double b =sqrt(pow(p.getX(), 2.0)+pow(p.getY(), 2.0)+pow(p.getZ(), 2.0));
        double angle = acos(dot_value/(a*b));
        return angle;
    }
};
Point pos; /// Declaration

class Ray{
    Point origin;
    Point direction;
public:
    Ray(){}
    Ray(Point origin, Point direction)
    {
        this->origin = origin;
        this->direction = direction;
        this->direction.normalize();
    }
    Point getOrigin(){return origin;}
    Point getDirection(){return direction;}
};
class PointLight{
    Point light_pos;
    Color clr;
    double radius;
    int stacks,slices;

public:
    PointLight(){
        light_pos.setPoint(0.0,0.0,0.0);
        clr.setColor(0.0,0.0,0.0);
        radius =0.0;stacks=0; slices=0;
    }
    PointLight(Point p, Color c, double radius, int slices, int stacks){
        this->light_pos=p; this->clr=c;this->radius = radius;this->slices=slices;this->stacks=stacks;
    }
    Point setPointAndColor(Point p, Color c){this->light_pos=p; this->clr=c;}
    Point getLightPosition(){return light_pos;}
    Color getColor(){return clr;}
    void drawOneSide(){
//         cout<<"into sphere"<<endl;
        Point points[stacks+1][slices+1];
        int i,j;
        double h,r;

        //generate points
        for(i=0;i<=stacks;i++){
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*(pi/2));
                points[i][j].y=r*sin(((double)j/(double)slices)*(pi/2));
                points[i][j].z=h;
            }
        }
        for(i=0;i<stacks;i++){
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                //upper hemisphere
				glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
				glVertex3f(points[i][j+1].x, points[i][j+1].y, points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x, points[i+1][j+1].y, points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x, points[i+1][j].y, points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
				glVertex3f(points[i][j+1].x, points[i][j+1].y, -points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x, points[i+1][j+1].y, -points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x, points[i+1][j].y, -points[i+1][j].z);
                }glEnd();
            }
        }
    }
    void drawSFourOfSphere(){
        ///upper hemisphere
        glPushMatrix();
        {
            /// +x,+y,+z
            drawOneSide();
        }
        glPopMatrix();
        glPushMatrix();
        {
            /// -x,+y,+z
            glRotated(90,0,0,1);
            drawOneSide();
        }
        glPopMatrix();
        glPushMatrix();
        {
            /// -x,-y,+z
            glRotated(180,0,0,1);
            drawOneSide();
        }
        glPopMatrix();
        glPushMatrix();
        {
            /// +x,-y,+z
            glRotated(270,0,0,1);
            drawOneSide();
        }
        glPopMatrix();

    }
    void draw(){

        glColor3f(clr.getR(), clr.getG(), clr.getB());
        glPushMatrix();
        {
            /// +x,+y,+z
            glTranslated(light_pos.getX(), light_pos.getY(), light_pos.getZ());
            drawSFourOfSphere();

        }
        glPopMatrix();
    }

};
class SpotLight{
    PointLight point_light;
    Point light_direction;
    double cutoff_angle;
    double radius;
    int stacks,slices;
public:
    SpotLight(){cutoff_angle=0.0;}
    SpotLight(PointLight pl, Point light_dir, double cutoff_angle, double radius, int slices, int stacks)
    {
        this->point_light = pl;
        this->light_direction= light_dir;
        this->cutoff_angle=cutoff_angle;
        this->radius = radius;
        this->slices=slices;
        this->stacks=stacks;
    }
    PointLight getPointLight(){return point_light;}
    Point getLightDirection(){return light_direction;}
    double getCuttoffAngle(){return cutoff_angle;}
    void drawOneSide(){
        Point points[stacks+1][slices+1];
        int i,j;
        double h,r;

        //generate points
        for(i=0;i<=stacks;i++){
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*(pi/2));
                points[i][j].y=r*sin(((double)j/(double)slices)*(pi/2));
                points[i][j].z=h;
            }
        }
        for(i=0;i<stacks;i++){
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                //upper hemisphere
				glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
				glVertex3f(points[i][j+1].x, points[i][j+1].y, points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x, points[i+1][j+1].y, points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x, points[i+1][j].y, points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
				glVertex3f(points[i][j+1].x, points[i][j+1].y, -points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x, points[i+1][j+1].y, -points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x, points[i+1][j].y, -points[i+1][j].z);
                }glEnd();
            }
        }
    }
    void drawSFourOfSphere(){
        ///upper hemisphere
        glPushMatrix();
        {
            /// +x,+y,+z
            drawOneSide();
        }
        glPopMatrix();
        glPushMatrix();
        {
            /// -x,+y,+z
            glRotated(90,0,0,1);
            drawOneSide();
        }
        glPopMatrix();
        glPushMatrix();
        {
            /// -x,-y,+z
            glRotated(180,0,0,1);
            drawOneSide();
        }
        glPopMatrix();
        glPushMatrix();
        {
            /// +x,-y,+z
            glRotated(270,0,0,1);
            drawOneSide();
        }
        glPopMatrix();

    }
    void draw(){


        glColor3f(point_light.getColor().getR(), point_light.getColor().getG(), point_light.getColor().getB());

        glPushMatrix();
        {
            /// +x,+y,+z
            glTranslated(point_light.getLightPosition().getX(), point_light.getLightPosition().getY(), point_light.getLightPosition().getZ());
            drawSFourOfSphere();

        }
        glPopMatrix();
    }
};
vector <PointLight> pointLights; ///Declaration
vector <SpotLight> spotLights;///Declaration

class Object{
protected:
    Point ref_point; // should have x, y, z
    double height, width, length;
    Color clr;
    double coEfficients[4]; /// ambient, diffuse,specular, reflection coefficients
    int shininess; /// exponent term of specular component
public:
    object(){
    }
    object(Point ref_point){this->ref_point = ref_point;}
    virtual void draw()=0;
    virtual double intersect(Ray,  int, Color&)=0;
    Color getColor() {return clr;}
    int getShininess(){return shininess;}
    void setColor(Color c){
        clr.setColor(c.getR(), c.getG(), c.getB());
    }
    void setShine(int s){ this->shininess = s;}
    void setCoEfficients(double a, double b, double c, double d){
        coEfficients[0]=a; ///ambient
        coEfficients[1]=b;///diffuse
        coEfficients[2]=c;///specular
        coEfficients[3]=d;///reflection coefficients
    }



    void ambientLight(Color& clr, Color int_point_clr) {
        double r =int_point_clr.getR()*coEfficients[0];
        double g = int_point_clr.getG()*coEfficients[0];
        double b = int_point_clr.getB()*coEfficients[0];
        clr.setColor(r,g,b);

    }
    void reflection(Ray ray, Color& clr, Point int_point, Color intersect_color, Point normal, PointLight point_light, Ray incidentRay) {
        Ray reflectedRay(int_point, incidentRay.getDirection()- normal*((normal.dotMultiplication(incidentRay.getDirection()))*2.0));
        double lambertValue =  normal.dotMultiplication(incidentRay.getDirection()*(-1.0));
        double phongValue = reflectedRay.getDirection().dotMultiplication((ray.getDirection()*(-1.0))); //Shine er jnno

        double max_val = max(lambertValue, 0.0);
        double r = point_light.getColor().getR()*intersect_color.getR()*coEfficients[1]*max_val;
        double g =point_light.getColor().getG()*intersect_color.getG()*coEfficients[1]*max_val;
        double b =point_light.getColor().getB()*intersect_color.getB()*coEfficients[1]*max_val;
        clr.setColor(clr.getR()+r, clr.getG()+g,clr.getB()+b);

        max_val = pow(max(phongValue, 0.0), getShininess());
        r =point_light.getColor().getR()*coEfficients[2]*max_val*intersect_color.getR() ;
        g =point_light.getColor().getG()*coEfficients[2]*max_val*intersect_color.getG() ;
        b =point_light.getColor().getB()*coEfficients[2]*max_val*intersect_color.getB() ;
        clr.setColor(clr.getR()+r, clr.getG()+g,clr.getB()+b);
    }
    void recursiveReflection(Color& clr, Color reflectedColor) {
        double r =reflectedColor.getR()*coEfficients[3];
        double g =reflectedColor.getG()*coEfficients[3];
        double b =reflectedColor.getB()*coEfficients[3];
        clr.setColor(clr.getR()+r, clr.getG()+g,clr.getB()+b);
    }
};
vector<Object*> objects; ///Declaration

class Sphere : public Object{
    Point center;
    double radius;
    int slices;
    int stacks;
//    Color c;
public:
    Sphere(Point center, double r, int slices,int stacks){

        this->center = center;
        this->radius = r;
        this->slices = slices;this->stacks = stacks;
    }
    Point getCenter(){return center;}
    Point setCenter(Point p){this->center = p;}
    void drawOneSide(){

        Point points[stacks+1][slices+1];
        int i,j;
        double h,r;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));

            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*(pi/2));
                points[i][j].y=r*sin(((double)j/(double)slices)*(pi/2));
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
            //glColor3f(100,0,0);
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
				glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
				glVertex3f(points[i][j+1].x, points[i][j+1].y, points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x, points[i+1][j+1].y, points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x, points[i+1][j].y, points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
				glVertex3f(points[i][j+1].x, points[i][j+1].y, -points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x, points[i+1][j+1].y, -points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x, points[i+1][j].y, -points[i+1][j].z);
                }glEnd();
            }
        }
    }
    void drawSFourOfSphere() {

        glPushMatrix();{
            drawOneSide();
        }
        glPopMatrix();
        glPushMatrix();{
            glRotated(90,0,0,1);
            drawOneSide();
        }
        glPopMatrix();
        glPushMatrix();{
            glRotated(180,0,0,1);
            drawOneSide();
        }
        glPopMatrix();
        glPushMatrix();{
            glRotated(270,0,0,1);
            drawOneSide();
        }
        glPopMatrix();
    }
    void draw(){
        glColor3f(getColor().getR(), getColor().getG(), getColor().getB());
        glPushMatrix();{
            glTranslated(center.getX(), center.getY(), center.getZ());
            drawSFourOfSphere();
        }
        glPopMatrix();
    }
    double intersect(Ray ray, int level, Color& clr){
        double a, b, c, tMin;

        double dir_dir, dir_org, dir_center, org_org, org_center, center_center;
        dir_dir = ray.getDirection().dotMultiplication(ray.getDirection()); //dir.dir
        dir_org = ray.getDirection().dotMultiplication(ray.getOrigin());
        dir_center = ray.getDirection().dotMultiplication(center);
        org_org = ray.getOrigin().dotMultiplication(ray.getOrigin());
        org_center = ray.getOrigin().dotMultiplication(center);
        center_center = center.dotMultiplication(center);
        ///Rd톀d t2 + 2Rd톀o t + Ro톀o - r2  =  0; a =Rd톀d , b=2Rd톀o, c= Ro톀o -r*r
        a = dir_dir;
        b = (dir_org-dir_center)*2.0;
        c = (org_org)+(center_center)-(org_center)*2.0- radius*radius;//

        double sqrt_val = b*b-4.0*a*c;

        if(sqrt_val < 0.0) {
            tMin = INF;
        }
        else if(sqrt_val > 0.0) {
            double tMax =( -b+sqrt(sqrt_val))/(2.0*a);
            tMin = ( -b-sqrt(sqrt_val))/(2.0*a);
            if(tMin < 0.0) tMin=tMax;///closest positive ta nite hbe
        }
        else {
            tMin = -b/(2.0*a);
        }

        if(level == 0) {
            return tMin;
        }
        Point int_Point = ray.getOrigin()+ray.getDirection()*tMin;
        Color intersect_color = getColor();

        Point normal = int_Point - center;
        normal.normalize();
        if(pos.distance(center) < radius) normal = normal*(-1.0);

        ambientLight(clr, intersect_color);
        ///Point light
        for(int i=0; i<pointLights.size(); i++) {
            Point origin = pointLights[i].getLightPosition();
            Point direction = int_Point - pointLights[i].getLightPosition();
            Ray inRay(origin, direction);
            double t, tMinimum=INF;
            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;
                t = objects[j]->intersect(inRay, 0, dummyColor);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }
            Point shadowPoint = inRay.getOrigin()+inRay.getDirection()*tMinimum;
            double dis1 = int_Point.distance(inRay.getOrigin());
            double dis2 = shadowPoint.distance(inRay.getOrigin());

            double epsilon = 0.0000001;
            if( (dis1 - epsilon) <  dis2)
            {
                reflection(ray, clr, int_Point, intersect_color, normal, pointLights[i], inRay);

            }
        }
        ///Spot lights
        for(int i=0; i<spotLights.size(); i++) {
            Point origin = spotLights[i].getPointLight().getLightPosition();
            Point direction = int_Point - spotLights[i].getPointLight().getLightPosition();
            Ray inRay(origin, direction);
            double t, tMinimum=INF;
            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;
                t = objects[j]->intersect(inRay, 0, dummyColor);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }
            double angle = spotLights[i].getLightDirection().getAngle(inRay.getDirection());
            double limit_angle =  spotLights[i].getCuttoffAngle()*(pi/180.0);

            if( angle <  limit_angle)
            {
                reflection(ray, clr, int_Point, intersect_color, normal, spotLights[i].getPointLight(), inRay);

            }
        }
        ///recursive
        if(level >= recur_level) {
            return tMin;
        }
        Point reflected_dir = ray.getDirection()- normal*((ray.getDirection().dotMultiplication(normal))*2.0);
        reflected_dir.normalize();
        Ray reflectedRay(int_Point+reflected_dir, reflected_dir);

        int near_idx = INT_MAX;
        double t, tempMin=INF;

        for(int i=0; i<objects.size(); i++) {
            Color dummyColor;
            t = objects[i]->intersect(reflectedRay, 0, dummyColor);

            if(t>0.0 && t<tempMin) {
                tempMin = t;
                near_idx = i;
            }
        }
        Color reflectedColor;

        if(near_idx != INT_MAX) {
            tempMin = objects[near_idx]->intersect(reflectedRay, level+1, reflectedColor);
        }
        recursiveReflection(clr, reflectedColor);

        ///color tuning
        double r,g,blue;
        if(clr.getR() > 1.0) r = 1.0;
        else
        {
            if(clr.getR() < 0.0) r = 0.0;
            else r = clr.getR();
        }

        if(clr.getG() > 1.0) g = 1.0;
       else
        {
            if(clr.getG() < 0.0) g = 0.0;
            else g = clr.getG();
        }

        if(clr.getB() > 1.0) blue = 1.0;
        else
        {
            if(clr.getB() < 0.0) blue = 0.0;
            else blue = clr.getB();
        }

        clr.setColor(r,g,blue);
        return tMin;
    }
};
class Triangle: public Object{
    Point a;
    Point b;
    Point c;
public:
    Triangle(Point a, Point b, Point c){
        this->a=a;
        this->b=b;
        this->c=c;
    }
    void draw()
    {
        glColor3f(getColor().getR(), getColor().getG(), getColor().getB());
         glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.getX(), a.getY(), a.getZ());
            glVertex3f(b.getX(), b.getY(), b.getZ());
            glVertex3f(c.getX(), c.getY(), c.getZ());
        }
        glEnd();
    }
    double intersect(Ray ray, int level, Color& color)
    {
        double A, beta, gamma, T;
        double dirx,diry,dirz,orgx,orgy,orgz;
        dirx = ray.getDirection().getX();
        diry = ray.getDirection().getY();
        dirz = ray.getDirection().getZ();
        orgx = ray.getOrigin().getX();
        orgy = ray.getOrigin().getY();
        orgz = ray.getOrigin().getZ();

        double first_col = (a.getY()-c.getY())*ray.getDirection().getZ()-(a.getZ()-c.getZ())*ray.getDirection().getY();
        double scnd_col = (a.getZ()-b.getZ())*ray.getDirection().getY()-(a.getY()-b.getY())*ray.getDirection().getZ();
        double third_col = (a.getY()-b.getY())*(a.getZ()-c.getZ())-(a.getZ()-b.getZ())*(a.getY()-c.getY());
        A = (a.getX()-b.getX())*(first_col) + (a.getX()-c.getX())*(scnd_col) + ray.getDirection().getX()*(third_col);

        first_col = (a.getY()-c.getY())*dirz-(a.getZ()-c.getZ())*diry;
        scnd_col = (a.getZ()-orgz)*diry-(a.getY()-orgy)*dirz;
        third_col = (a.getY()-orgy)*(a.getZ()-c.getZ())-(a.getZ()-orgz)*(a.getY()-c.getY());

        beta = (a.getX()-orgx)*(first_col) + (a.getX()-c.getX())*(scnd_col) + dirx*(third_col);

        first_col = (a.getY()-ray.getOrigin().getY())*ray.getDirection().getZ()-(a.getZ()-ray.getOrigin().getZ())*ray.getDirection().getY();
        scnd_col = (a.getZ()-b.getZ())*ray.getDirection().getY()-(a.getY()-b.getY())*ray.getDirection().getZ();
        third_col = (a.getY()-b.getY())*(a.getZ()-ray.getOrigin().getZ())-(a.getZ()-b.getZ())*(a.getY()-ray.getOrigin().getY());

        gamma = (a.getX()-b.getX())*(first_col) + (a.getX()-ray.getOrigin().getX())*(scnd_col) + ray.getDirection().getX()*(third_col);

        first_col = (a.getY()-c.getY())*(a.getZ()-ray.getOrigin().getZ())-(a.getZ()-c.getZ())*(a.getY()-ray.getOrigin().getY());
        scnd_col = (a.getZ()-b.getZ())*(a.getY()-ray.getOrigin().getY())-(a.getY()-b.getY())*(a.getZ()-ray.getOrigin().getZ());
        third_col = (a.getY()-b.getY())*(a.getZ()-c.getZ())-(a.getZ()-b.getZ())*(a.getY()-c.getY());

        T = (a.getX()-b.getX())*(first_col);
        T += (a.getX()-c.getX())*(scnd_col);
        T += (a.getX()-ray.getOrigin().getX())*(third_col);

        double tMin=INF;
        if(A != 0.0) {

            if(beta/A>0.0 && gamma/A>0.0 && beta/A+gamma/A<1.0) {
                tMin = T/A;
            }
            else {
                tMin = INF;
            }
        }
        if(level == 0) {
            return tMin;
        }
        Point int_Point = ray.getOrigin()+ray.getDirection()*tMin;
        Color intersect_color = getColor();

        Point m = b-a;
        Point n=c-a;
        Point normal = m.crossMultiplication(n);
        normal.normalize();
        normal = ((normal.dotMultiplication(ray.getDirection()*(-1.0))) > 0.0)? normal: normal*(-1.0);

        ambientLight(color, intersect_color);
        ///Point light
        for(int i=0; i<pointLights.size(); i++) {
            Point origin = pointLights[i].getLightPosition();
            Point direction = int_Point - pointLights[i].getLightPosition();
            Ray inRay(origin, direction);
            double t, tMinimum=INF;
            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;
                t = objects[j]->intersect(inRay, 0, dummyColor);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }
            Point shadowPoint = inRay.getOrigin()+inRay.getDirection()*tMinimum;
            double dis1 = int_Point.distance(inRay.getOrigin());
            double dis2 = shadowPoint.distance(inRay.getOrigin());

            double epsilon = 0.0000001;
            if( (dis1 - epsilon) <  dis2)
            {
                reflection(ray, color, int_Point, intersect_color, normal, pointLights[i], inRay);

            }
        }
        ///Spot lights
        for(int i=0; i<spotLights.size(); i++) {
            Point origin = spotLights[i].getPointLight().getLightPosition();
            Point direction = int_Point - spotLights[i].getPointLight().getLightPosition();
            Ray inRay(origin, direction);
            double t, tMinimum=INF;
            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;
                t = objects[j]->intersect(inRay, 0, dummyColor);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }
            double angle = spotLights[i].getLightDirection().getAngle(inRay.getDirection());
            double limit_angle =  spotLights[i].getCuttoffAngle()*(pi/180.0);

            if( angle <  limit_angle)
            {
                reflection(ray, color, int_Point, intersect_color, normal, spotLights[i].getPointLight(), inRay);

            }
        }
        ///recursive
        if(level >= recur_level) {
            return tMin;
        }

        double temp = ray.getDirection().dotMultiplication(normal);
        Point reflected_dir = ray.getDirection()- normal*temp*2.0;
        reflected_dir.normalize();
        Ray reflectedRay(int_Point+reflected_dir, reflected_dir);

        int near_idx = INT_MAX;
        double t, tempMin=INF;

        for(int i=0; i<objects.size(); i++) {
            Color dummyColor;
            t = objects[i]->intersect(reflectedRay, 0, dummyColor);

            if(t>0.0 && t<tempMin) {
                tempMin = t;
                near_idx = i;
            }
        }
        Color reflectedColor;

        if(near_idx != INT_MAX) {
            tempMin = objects[near_idx]->intersect(reflectedRay, level+1, reflectedColor);
        }
        recursiveReflection(color, reflectedColor);

        double r,g,blue;
        if(color.getR() > 1.0) r = 1.0;

        else
        {
            if(color.getR() < 0.0) r = 0.0;
            else r = color.getR();
        }

        if(color.getG() > 1.0) g = 1.0;
       else
        {
            if(color.getG() < 0.0) g = 0.0;
            else g = color.getG();
        }

        if(color.getB() > 1.0) blue = 1.0;
        else
        {
            if(color.getB() < 0.0) blue = 0.0;
            else blue = color.getB();
        }

        color.setColor(r,g,blue);

        return tMin;
    }

};
class General: public Object{
    GeneralPoints coeff;
    Point cube_ref_point;
    double length,width,height;
public:
    General() {
        length = width = height = 0.0;
    }
    General(struct GeneralPoints coeff, Point cube_ref_point, double length, double width, double height) {
        this->coeff = coeff;
        this->cube_ref_point = cube_ref_point;
        this->length= length;
        this->width = width;
        this->height = height;
    }

    void draw() {}
    double intersect(Ray ray, int level, Color& clr)
    {
        double a, b, c;
        double dirx_dirx, diry_diry, dirz_dirz, dirx_diry, dirx_dirz, diry_dirz;
        dirx_dirx = ray.getDirection().getX()*ray.getDirection().getX();
        diry_diry = ray.getDirection().getY()*ray.getDirection().getY();
        dirz_dirz = ray.getDirection().getZ()*ray.getDirection().getZ();
        dirx_diry = ray.getDirection().getX()*ray.getDirection().getY();
        dirx_dirz = ray.getDirection().getX()*ray.getDirection().getZ();
        diry_dirz = ray.getDirection().getY()*ray.getDirection().getZ();

        a = coeff.a*dirx_dirx + coeff.b*diry_diry +coeff.c*dirz_dirz + coeff.d*dirx_diry+ coeff.e*dirx_dirz+coeff.f*diry_dirz;

        double dirx,diry,dirz,orgx,orgy,orgz;
        dirx = ray.getDirection().getX();
        diry = ray.getDirection().getY();
        dirz = ray.getDirection().getZ();
        orgx = ray.getOrigin().getX();
        orgy = ray.getOrigin().getY();
        orgz = ray.getOrigin().getZ();

        b = 2.0*coeff.a*orgx*dirx+2.0*coeff.b*orgy*diry+2.0*coeff.c*orgz*dirz;
        b += coeff.d*(orgx*diry+dirx*orgy);
        b += coeff.e*(orgx*dirz+dirx*orgz);
        b += coeff.f*(orgy*dirz+diry*orgz);
        b += coeff.g*dirx + coeff.h*diry + coeff.i*dirz;

        c = coeff.a*orgx*orgx + coeff.b*orgy*orgy+coeff.c*orgz*orgz ;
        c += coeff.d*orgx*orgy+coeff.e*orgx*orgz+coeff.f*orgy*orgz;
        c += coeff.g*orgx+coeff.h*orgy+coeff.i*orgz+coeff.j;


        double  tMin = INF, tMax =INF;
        double sqrt_val = b*b-4.0*a*c;

        if(a == 0.0)
        {
            if(b == 0.0) tMin = INF;
            else tMin=-c/b;
            tMax = INF;
        }
        else
        {

            if(sqrt_val < 0.0) {
                tMin = INF;
                tMax = INF;
            }
            else if(sqrt_val > 0.0)
            {
                tMax = (-b+sqrt(sqrt_val))/(2.0*a);
                tMin = (-b-sqrt(sqrt_val))/(2.0*a);
//                if(tMin < 0.0) tMin=tMax;
            }
            else {
                tMin = -b/(2.0*a);
                tMax = INF;
            }
        }
        if(tMin < INF)
        {
            if(tMax < INF)
            {

                if(tMin > 0.0)
                {
                    Point int_Point = ray.getOrigin()+ray.getDirection()*tMin;
                    int len_check = length!=0.0 &&(int_Point.getX()<cube_ref_point.getX() || (int_Point.getX()>cube_ref_point.getX()+length));
                    int width_check = width!=0.0 && (int_Point.getY()<cube_ref_point.getY() || (int_Point.getY()>cube_ref_point.getY()+width));
                    int height_check = height!=0.0 && (int_Point.getZ()<cube_ref_point.getZ() || (int_Point.getZ()>cube_ref_point.getZ()+height));
                    if(len_check)
                    {
                        tMin = INF;
                    }
                    else if(width_check)
                    {
                        tMin = INF;
                    }
                    else if (height_check)
                    {
                        tMin = INF;
                    }

                }
                if(tMax > 0.0)
                {
                    Point int_Point = ray.getOrigin()+ray.getDirection()*tMax;
                    int len_check = length!=0.0 &&(int_Point.getX()<cube_ref_point.getX() || (int_Point.getX()>cube_ref_point.getX()+length));
                    int width_check = width!=0.0 && (int_Point.getY()<cube_ref_point.getY() || (int_Point.getY()>cube_ref_point.getY()+width));
                    int height_check = height!=0.0 && (int_Point.getZ()<cube_ref_point.getZ() || (int_Point.getZ()>cube_ref_point.getZ()+height));

                   if(len_check)
                    {
                        tMax = INF;
                    }
                    else if(width_check)
                    {
                        tMax = INF;
                    }
                    else if (height_check)
                    {
                        tMax = INF;
                    }
                }
                if(tMin>0.0 && tMin<tMax) tMin = tMin;
                else tMin =tMax;
            }
           else
            {
                if(tMin > 0.0) {
                    Point int_Point = ray.getOrigin()+ray.getDirection()*tMin;
                    int len_check = length!=0.0 &&(int_Point.getX()<cube_ref_point.getX() || (int_Point.getX()>cube_ref_point.getX()+length));
                    int width_check = width!=0.0 && (int_Point.getY()<cube_ref_point.getY() || (int_Point.getY()>cube_ref_point.getY()+width));
                    int height_check = height!=0.0 && (int_Point.getZ()<cube_ref_point.getZ() || (int_Point.getZ()>cube_ref_point.getZ()+height));
                    if(len_check)
                    {
                        tMin = INF;
                    }
                    else if(width_check)
                    {
                        tMin = INF;
                    }
                    else if (height_check)
                    {
                        tMin = INF;
                    }
                }
            }
        }

        if(level == 0) {return tMin;}
        double xN, yN, zN;
        Point int_Point = ray.getOrigin()+ray.getDirection()*tMin;
        Color intersect_color = getColor();

        xN = 2.0*coeff.a*int_Point.getX()+coeff.d*int_Point.getY() + coeff.e*int_Point.getZ()+coeff.g;
        yN = 2.0*coeff.b*int_Point.getY()+coeff.d*int_Point.getX()+coeff.f*int_Point.getZ()+coeff.h;
        zN = 2.0*coeff.c*int_Point.getZ()+coeff.e*int_Point.getX() + coeff.f*int_Point.getY()+coeff.i;

        Point normal(xN, yN, zN);
        normal.normalize();
        if(normal.dotMultiplication(ray.getDirection()*(-1.0)) <0.0)normal = normal*(-1.0);
        ambientLight(clr, intersect_color);

        for(int i=0; i<pointLights.size(); i++) {
            Point origin = pointLights[i].getLightPosition();
            Point direction = int_Point - pointLights[i].getLightPosition();
            Ray inRay(origin, direction);
            double t, tMinimum=INF;
            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;
                t = objects[j]->intersect(inRay, 0, dummyColor);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }
            Point shadowPoint = inRay.getOrigin()+inRay.getDirection()*tMinimum;
            double dis1 = int_Point.distance(inRay.getOrigin());
            double dis2 = shadowPoint.distance(inRay.getOrigin());

            double epsilon = 0.0000001;
            if( (dis1 - epsilon) <  dis2)
            {
                reflection(ray, clr, int_Point, intersect_color, normal, pointLights[i], inRay);

            }
        }
        ///Spot lights
        for(int i=0; i<spotLights.size(); i++) {
            Point origin = spotLights[i].getPointLight().getLightPosition();
            Point direction = int_Point - spotLights[i].getPointLight().getLightPosition();
            Ray inRay(origin, direction);
            double t, tMinimum=INF;
            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;
                t = objects[j]->intersect(inRay, 0, dummyColor);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }
            double angle = spotLights[i].getLightDirection().getAngle(inRay.getDirection());
            double limit_angle =  spotLights[i].getCuttoffAngle()*(pi/180.0);

            if( angle <=  limit_angle)
            {
                reflection(ray, clr, int_Point, intersect_color, normal, spotLights[i].getPointLight(), inRay);

            }
        }
        ///recursive
        if(level >= recur_level) {
            return tMin;
        }

        double temp = ray.getDirection().dotMultiplication(normal);
        Point reflected_dir = ray.getDirection()- normal*temp*2.0;
        reflected_dir.normalize();
        Ray reflectedRay(int_Point+reflected_dir, reflected_dir);

        int near_idx = INT_MAX;
        double t, tempMin=INF;

        for(int i=0; i<objects.size(); i++) {
            Color dummyColor;
            t = objects[i]->intersect(reflectedRay, 0, dummyColor);

            if(t>0.0 && t<tempMin) {
                tempMin = t;
                near_idx = i;
            }
        }
        Color reflectedColor;

        if(near_idx != INT_MAX) {
            tempMin = objects[near_idx]->intersect(reflectedRay, level+1, reflectedColor);
        }
        recursiveReflection(clr, reflectedColor);

         double red,green,blue;
        if(clr.getR() > 1.0) red = 1.0;

        else
        {
            if(clr.getR() < 0.0) red = 0.0;
            else red = clr.getR();
        }

        if(clr.getG() > 1.0) green = 1.0;
       else
        {
            if(clr.getG() < 0.0) green = 0.0;
            else green = clr.getG();
        }

        if(clr.getB() > 1.0) blue = 1.0;
        else
        {
            if(clr.getB() < 0.0) blue = 0.0;
            else blue = clr.getB();
        }

        clr.setColor(red,green,blue);
//
        return tMin;
    }
};
class Floor: public Object{
    int floorWidth, tileWidth;
    Color clr;
public:
    Floor() {
        floorWidth = tileWidth = 0.0;
    }
    Floor(int floorWidth, int tileWidth, Color clr){
        ref_point.setX(-floorWidth/2);
        ref_point.setY(-floorWidth/2);
        ref_point.setZ(0.0);
        this->tileWidth=tileWidth;
        this->floorWidth=floorWidth;
        this->clr = clr;
    }
    void draw(){
        double r=0.0, g =0, b=0;
        for(int i=0; i< (floorWidth/tileWidth); i++){

            for(int j=0; j< (floorWidth/tileWidth); j++){
                if((i+j)%2 == 0)
                {
                    glColor3f(getColor().getR(),getColor().getG(), getColor().getB());
                }
                else
                {
                    glColor3f(clr.getR(),clr.getG(), clr.getB());
                }
                double x = -floorWidth/2.0+tileWidth*j ;
                double y = -floorWidth/2.0+tileWidth*i ;
                double z = 0.0;
                glBegin(GL_QUADS);
                {
                    glVertex3f(x, y, z);
                    glVertex3f(x+tileWidth, y, z);
                    glVertex3f(x+tileWidth, y+tileWidth, z);
                    glVertex3f(x, y+tileWidth, z);
                }
                glEnd();
;

            }

        }
    }
    double intersect(Ray ray, int level, Color& color){
        Point normal(0.0, 0.0, 1.0);
        double check = normal.dotMultiplication(pos);
        if(check < 0.0) normal= normal*(-1.0); // - z
        double tMin = INF;
        double tCurrent;

        check = normal.dotMultiplication(ray.getDirection()); // check whether it is divided by 0 or not

        if( normal.dotMultiplication(ray.getDirection()) != 0.0) {
            double a = normal.dotMultiplication(ray.getOrigin());
            double b = normal.dotMultiplication(ray.getDirection());
            tMin = (-1.0)*(a/b);
        }

        double floorLen = floorWidth/2.0;

        if(tMin>0.0 && tMin<INF) {
            Point int_point = ray.getOrigin() + ray.getDirection()*tMin;
                if(!((int_point.getX()>-floorLen && int_point.getX()<floorLen ) && (int_point.getY()>-floorLen  && int_point.getY()<floorLen ))) {
                    tMin = INF;
                }
        }


        if(level == 0) {

            return tMin;
        }
        Point int_Point = ray.getOrigin() + ray.getDirection()*tMin;
        Point referencePosition = int_Point - Point(-floorLen, -floorLen, 0.0);
        Color intersect_color = (((int) (floor(referencePosition.getX()/tileWidth)+
                                                floor(referencePosition.getY()/tileWidth)))%2 == 0)
                                                ? getColor()
                                                : clr;

        ambientLight(color, intersect_color);
        ///Point light
        for(int i=0; i<pointLights.size(); i++) {
            Point origin = pointLights[i].getLightPosition();
            Point direction = int_Point - pointLights[i].getLightPosition();
            Ray inRay(origin, direction);
            double t, tMinimum=INF;
            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;
                t = objects[j]->intersect(inRay, 0, dummyColor);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }
            Point shadowPoint = inRay.getOrigin()+inRay.getDirection()*tMinimum;
            double dis1 = int_Point.distance(inRay.getOrigin());
            double dis2 = shadowPoint.distance(inRay.getOrigin());

            double epsilon = 0.0000001;
            if( (dis1 - epsilon) <  dis2)
            {
                reflection(ray, color, int_Point, intersect_color, normal, pointLights[i], inRay);

            }
        }
        ///Spot lights
        for(int i=0; i<spotLights.size(); i++) {
            Point origin = spotLights[i].getPointLight().getLightPosition();
            Point direction = int_Point - spotLights[i].getPointLight().getLightPosition();
            Ray inRay(origin, direction);
            double t, tMinimum=INF;
            for(int j=0; j<objects.size(); j++) {
                Color dummyColor;
                t = objects[j]->intersect(inRay, 0, dummyColor);

                if(t>0.0 && t<tMinimum) {
                    tMinimum = t;
                }
            }
            double angle = spotLights[i].getLightDirection().getAngle(inRay.getDirection());
            double limit_angle =  spotLights[i].getCuttoffAngle()*(pi/180.0);

            if( angle <  limit_angle)
            {
                reflection(ray, color, int_Point, intersect_color, normal, spotLights[i].getPointLight(), inRay);

            }
        }
        ///recursive
        if(level >= recur_level) {
            return tMin;
        }
        Point reflected_dir = ray.getDirection()- normal*((ray.getDirection().dotMultiplication(normal))*2.0);
        reflected_dir.normalize();
        Ray reflectedRay(int_Point+reflected_dir, reflected_dir);

        int near_idx = INT_MAX;
        double t, tempMin=INF;

        for(int i=0; i<objects.size(); i++) {
            Color dummyColor;
            t = objects[i]->intersect(reflectedRay, 0, dummyColor);

            if(t>0.0 && t<tempMin) {
                tempMin = t;
                near_idx = i;
            }
        }
        Color reflectedColor;

        if(near_idx != INT_MAX) {
            tempMin = objects[near_idx]->intersect(reflectedRay, level+1, reflectedColor);
        }
        recursiveReflection(color, reflectedColor);

        ///color tuning
        double r,g,blue;
        if(color.getR() > 1.0) r = 1.0;
        else
        {
            if(color.getR() < 0.0) r = 0.0;
            else r = color.getR();
        }

        if(color.getG() > 1.0) g = 1.0;
       else
        {
            if(color.getG() < 0.0) g = 0.0;
            else g = color.getG();
        }

        if(color.getB() > 1.0) blue = 1.0;
        else
        {
            if(color.getB() < 0.0) blue = 0.0;
            else blue = color.getB();
        }

        color.setColor(r,g,blue);

        return tMin;
    }

};





