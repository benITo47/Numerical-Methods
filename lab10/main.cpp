#include <iostream> 
#include <array>
#include <cmath> 

double func(double x, double y){

    double result = sin(x) * sin(y) - exp(-((x+M_PI/2)*(x+M_PI/2)) - ((y - M_PI/2)* (y - M_PI/2)));

    return result; 
}

double d_rand ( const double min , const double max ){
 double r = ( double ) rand () / RAND_MAX ; // Przedzial [0 , 1]
 r = r * ( max - min ) + min ; // Przeskalowanie do [min , max]
 return r;
}

double checkBounds(double val, double min, double max)
{
    if(val > max)
    {
        val = max;
    }
    else if(val < min)
    {
        val = min; 
    }
    return val;
}

int main()
{
    std::array<std::array<double,2>,200> wedrowcy;

    FILE* plik1 = fopen("w0.dat", "w");
    FILE* plik2 = fopen("T.dat", "w");
    

    double x_min, y_min;
    double x_max, y_max;
    x_min = y_min = -10.0; 
    x_max = y_max = 10.0; 

    double x_start, y_start; 
    x_start = y_start = 5;

    for(int i = 0; i < 200; i++)
    {
        wedrowcy[i][0] = x_start; 
        wedrowcy[i][1] = y_start; 
    }

        for(int temp = 0; temp <= 20; temp++)
        {
            double Temperatura = 10/(pow(2,temp));
            for(int k = 0; k < 100; k++)
            {
                for(int w = 0; w < 200; w++){
                double x = wedrowcy[w][0];
                double y = wedrowcy[w][1];
                double delta_x = d_rand(-1,1);
                double delta_y = d_rand(-1,1);

                double x_dx = x + delta_x; 
                double y_dy = y + delta_y; 

                x_dx = checkBounds(x_dx,x_min,x_max);
                y_dy = checkBounds(y_dy,y_min,y_max);

                if(func(x_dx,y_dy) < func(x,y))
                {
                    wedrowcy[w][0] = x_dx; 
                    wedrowcy[w][1] = y_dy;
                }
                else if(d_rand(0,1) < exp(- (func(x_dx, y_dy) - func(x,y))/Temperatura))
                {
                    wedrowcy[w][0] = x_dx; 
                    wedrowcy[w][1] = y_dy;
                }
                }

                fprintf(plik1, "%g  \n", func(wedrowcy[0][0], wedrowcy[0][1]) );
            }

           
            if(temp == 0 || temp == 7 || temp == 20)
                 {
                      for(int i = 0 ; i < 200; i++)
                         fprintf(plik2, "%g %g \n",wedrowcy[i][0], wedrowcy[i][1]);

                    fprintf(plik2, "\n\n\n\n\n\n\n\n\n\n\n" ); 

                }


        }
    
        //szukamy minimum, po wedrowcach

        double minimum = 1000;
        int i_min = 0; 
        for(int i = 0; i < 200; i++)
        {
            if(func(wedrowcy[i][0], wedrowcy[i][1]) < minimum)
            {
                i_min = i; 
                minimum = func(wedrowcy[i][0], wedrowcy[i][1]);
            }
        }

        std::cout << "\n\n\nMimimum znaleziono w x: " << wedrowcy[i_min][0] << " y:" << wedrowcy[i_min][1] << " i wynosi: " << minimum << std::endl;


    fclose(plik1);
    fclose(plik2);
}