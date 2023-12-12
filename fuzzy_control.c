#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <dir.h>
#include "libraries/pbPlots.c"
#include "libraries/supportLib.c"

#define abs_float(x) ((x)>0?(x):-(x))

// Global Variables
char rules[2][2][3];
char rules_list[5][3] = {"NB", "NS", "ZE", "PS", "PB"};
double m[2][2];
double w[4];
double force[4];

// Function to add char in a string
char add(char *s, const int len, char x, int pos){
    s[len] = 0;
    memmove(s+pos+1, s+pos, len-pos+1);
    s[pos]=x;
}

// Function to find the index from rules_list
int find_idx(char *x){
    int idx;
    for(idx = 0; idx<5; idx++){
        if(!strcmp(rules_list[idx], x)){
            break;
        }
    }

    return idx;
}

// Calculations for rise time
double rise_time(double time[], double y[], int iteration){
    double y_target = 0.1;
    double ten_percent;
    double ninety_percent;

    for(int i = 0; i<2; i++){
        int idx;
        for(idx=0; idx<iteration; idx++){
            if(y[idx] > y_target){
                break;
            }
        }

        double xtr1 = time[idx-1];
        double xtr2 = time[idx];
        double ytr1 = y[idx-1];
        double ytr2 = y[idx];
        double grad = (ytr2 - ytr1) / (xtr2 - xtr1);
        double y_intercept = ytr1-(grad*xtr1);

        if(i == 0)
            ten_percent = ((y_target-y_intercept)/grad);
        else
            ninety_percent = ((y_target-y_intercept)/grad);

        y_target = 0.9;
    }

    return (ninety_percent - ten_percent);
}

// Calculations for steady state error
double steady_state_error(double y[], double d[], int iteration){
    double sse = INFINITY;
    for(int i = 1; i<iteration; i++){
        if(abs_float(y[i] - y[i-1]) < 1e-5){
            sse = d[i] - y[i];
            break;
        }
    }

    return sse;
}

// Calculations for overshoot
double overshoot(double y[], double d[], int iteration){
    double os = 0;
    for(int i = 0; i<iteration; i++){
        if(y[i] > os && y[i] > d[i]){
            os = y[i];
        }
    }

    if(os != 0)
        os -= d[0];

    return os;
}

// Function to calculate mu
void calculate_mu(double error, double a, double b, int idx){
    double m1 = (b - error)/(b-a);
    double m2 = (error - a)/(b-a);

    m[idx][0] = m1;
    m[idx][1] = m2;
}

// Find rules and calculate mu
void membership_function(double error, double e[], int idx){
    double e1 = e[0];
    double e2 = e[1];
    double e3 = e[2];
    double e4 = e[3];
    double e5 = e[4];

    if(error <= e1){
        m[idx][0] = 1;
        m[idx][1] = 0;
        rules[idx][0][0] = 'N';
        rules[idx][0][1] = 'B';
        rules[idx][1][0] = 'N';
        rules[idx][1][1] = 'S';
    }else if(error < e2){
        calculate_mu(error, e1, e2, idx);
        rules[idx][0][0] = 'N';
        rules[idx][0][1] = 'B';
        rules[idx][1][0] = 'N';
        rules[idx][1][1] = 'S';
    }else if(error == e2){
        m[idx][0] = 0;
        m[idx][1] = 1;
        rules[idx][0][0] = 'N';
        rules[idx][0][1] = 'B';
        rules[idx][1][0] = 'N';
        rules[idx][1][1] = 'S';
    }else if(error < e3){
        calculate_mu(error, e2, e3, idx);
        rules[idx][0][0] = 'N';
        rules[idx][0][1] = 'S';
        rules[idx][1][0] = 'Z';
        rules[idx][1][1] = 'E';
    }else if(error == e3){
        m[idx][0] = 0;
        m[idx][1] = 1;
        rules[idx][0][0] = 'N';
        rules[idx][0][1] = 'S';
        rules[idx][1][0] = 'Z';
        rules[idx][1][1] = 'E';
    }else if(error < e4){
        calculate_mu(error, e3, e4, idx);
        rules[idx][0][0] = 'Z';
        rules[idx][0][1] = 'E';
        rules[idx][1][0] = 'P';
        rules[idx][1][1] = 'S';
    }else if(error == e4){
        m[idx][0] = 0;
        m[idx][1] = 1;
        rules[idx][0][0] = 'Z';
        rules[idx][0][1] = 'E';
        rules[idx][1][0] = 'P';
        rules[idx][1][1] = 'S';
    }else if(error < e5){
        calculate_mu(error, e4, e5, idx);
        rules[idx][0][0] = 'P';
        rules[idx][0][1] = 'S';
        rules[idx][1][0] = 'P';
        rules[idx][1][1] = 'B';
    }else if(error >= e5){
        m[idx][0] = 0;
        m[idx][1] = 1;
        rules[idx][0][0] = 'P';
        rules[idx][0][1] = 'S';
        rules[idx][1][0] = 'P';
        rules[idx][1][1] = 'B';
    }
}

void fuzzy_control(double defuzz[], double err_range[], double delta_err_range[], int num_plant, int iteration){
    char path[100] = "outputs/plant.png";
    char plant[2];
    char plant_str[10] = "PLANT ";
    wchar_t plant_str2[10];

    sprintf(plant, "%d", num_plant);

    add(path, 100, *plant, 13);
    add(plant_str, 10, *plant, 6);

    mbstowcs(plant_str2, plant_str, 10);

    printf("\n===========%s===========\n", plant_str);

    // Output Variable
    double y[500] = {0.0};
    double x2[500] = {0.0};
    double x1[500] = {0.0};
    double y_delta_u[500] = {};
    
    // Fuzzy output variable
    double u[500] = {0.0};
    // Desired Output
    double d[500] = {0.0};
    // Time variable
    double time[500];
    // Previous error
    double err_pref = 0.0;
    // Current error
    double err = d[0] - y[0];
    double delta_err = 0.0;

    // Bool to switch between plotting Δu(k) or u(k) 
    bool is_delta_uk = true;

    // Get all defuzzifier value
    double NB = defuzz[0];
    double NS = defuzz[1];
    double ZE = defuzz[2];
    double PS = defuzz[3];
    double PB = defuzz[4];

    // Force of rules matrix
    double rules_matrix[5][5] = {{NB, NB, NB, NS, PB},
                                 {NB, NS, NS, ZE, PB},
                                 {NB, NS, ZE, PS, PB},
                                 {NB, ZE, PS, PS, PB},
                                 {NB, PS, PB, PB, PB}};

    // Set value for time and desired output
    for(int i = 0; i<iteration; i++){
        time[i] = (double)i;
        d[i] = 1;
    }

    for(int z = 0; z<2; z++){
        for(int k = 0; k<iteration-1; k++){
            // Error calculations
            err_pref = err;
            err = d[k] - y[k];
            delta_err = err - err_pref;

            // Fuzzifier
            membership_function(err, err_range, 0);
            membership_function(delta_err, delta_err_range, 1);

            // Calculate weights
            for(int i = 0; i<2; i++){
                for(int j = 0; j<2; j++){
                    w[(i*2)+j] = m[0][i]*m[1][j];
                }
            }

            // Defuzzifier
            for(int i = 0; i<2; i++){
                int i_idx = find_idx(rules[0][i]);
                
                for(int j = 0; j<2; j++){
                    int j_idx = find_idx(rules[1][j]);

                    force[(i*2)+j] = rules_matrix[j_idx][i_idx];
                }
            }

            u[k+1] = 0;
            double w_tot = 0;
            
            // Calculate fuzzy output
            for(int i = 0; i<4; i++){
                w_tot += w[i];
                u[k+1] += w[i]*force[i];
            }
            u[k+1] /= w_tot;

            // If fuzzy output as Δu(k)
            if(is_delta_uk){
                u[k+1] = u[k] + u[k+1];
            }

            if(num_plant == 1){
                // 1st Plant
                y[k+1] = ((1.01*y[k]) + (0.01*(pow(y[k], 2))) + (0.03*u[k+1]));
            }else if(num_plant == 2){
                // 2nd Plant
                x1[k+1] = (x1[k]+(0.01*x2[k])+(0.01*u[k+1]));
                x2[k+1] = ((0.1*x1[k])+(0.97*x2[k]));
                y[k+1] = x1[k+1];
            }
        }

        // Calculate Rise Time
        double tr = rise_time(time, y, iteration);

        // Calculate SSE
        double sse = steady_state_error(y, d, iteration);

        // Calculate Overshoot
        double os = overshoot(y, d, iteration);

        if(is_delta_uk)
            printf("%s", "\nFuzzy Output as delta_u(k):");
        else
            printf("%s", "\nFuzzy Output as u(k):");
        printf("\n Rise Time : %f", tr);
        printf("\n Overshoot : %f", os);
        if(sse != INFINITY)
            printf("\n SSE       : %f\n", sse);
        else
            printf("\n SSE       : %s", "System is Unsteady\n");

        // Transition to fuzzy output as u(k) 
        if(is_delta_uk){
            // Set variable y_Δu(k)
            for(int i = 0; i<iteration; i++){
                y_delta_u[i] = y[i];
            }

            is_delta_uk = false;
        }
    }

    // Create directory for saving plot results
    mkdir("outputs");

    // Plot when fuzzy output is u(k) 
    ScatterPlotSeries *series = GetDefaultScatterPlotSeriesSettings();
	series->xs = time;
	series->xsLength = iteration;
	series->ys = y;
	series->ysLength = iteration;
	series->linearInterpolation = true;
	series->lineType = L"solid";
	series->lineTypeLength = wcslen(series->lineType);
	series->lineThickness = 2;
	series->color = CreateRGBColor(0, 0, 1);

    // Plot when fuzzy output is Δu(k)
    ScatterPlotSeries *series2 = GetDefaultScatterPlotSeriesSettings();
	series2->xs = time;
	series2->xsLength = iteration;
	series2->ys = y_delta_u;
	series2->ysLength = iteration;
	series2->linearInterpolation = true;
	series2->lineType = L"solid";
	series2->lineTypeLength = wcslen(series2->lineType);
	series2->lineThickness = 2;
	series2->color = CreateRGBColor(1, 0, 0);

    // Plot desired output
    ScatterPlotSeries *series3 = GetDefaultScatterPlotSeriesSettings();
	series3->xs = time;
	series3->xsLength = iteration;
	series3->ys = d;
	series3->ysLength = iteration;
	series3->linearInterpolation = true;
	series3->lineType = L"dashed";
	series3->lineTypeLength = wcslen(series3->lineType);
	series3->lineThickness = 2;
	series3->color = GetGray(0.6);

    // Plot settings
	ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
	settings->width = 1000;
	settings->height = 600;
	settings->autoBoundaries = true;
	settings->autoPadding = true;
    settings->title = plant_str2;
    settings->titleLength = wcslen(settings->title);
    settings->xLabel = L"Y(k)";
    settings->xLabelLength = wcslen(settings->xLabel);
    settings->yLabel = L"Steps";
    settings->yLabelLength = wcslen(settings->yLabel);
	ScatterPlotSeries *s [] = {series3, series2, series};
	settings->scatterPlotSeries = s;
	settings->scatterPlotSeriesLength = 3;

	RGBABitmapImageReference *canvasReference = CreateRGBABitmapImageReference();
	DrawScatterPlotFromSettings(canvasReference, settings);

	size_t length;
	double *pngdata = ConvertToPNG(&length, canvasReference->image);
	WriteToFile(pngdata, length, path);
	DeleteImage(canvasReference->image);

    printf("\n===============================\n", plant_str);
    printf("\n[INFO] Plot Saved To: '%s'\n", path);
}

int main(){
    // Defuzzifier value
    double uk[5] = {-10, -5, 0, 5, 10};

    // Error initial parameter value
    double ek[5] = {-1, -0.5, 0, 0.5, 1};

    // Delta error initial parameter value
    double delta_ek[5] = {-1, -0.5, 0, 0.5, 1};

    // Calculate and Plot 1st Plant
    delta_ek[0] = -0.35;
    delta_ek[1] = -0.14;
    delta_ek[2] = 0;
    delta_ek[3] = 0.14;
    delta_ek[4] = 0.35;
    fuzzy_control(uk, ek, delta_ek, 1, 100);

    // Calculate and Plot 2nd Plant
    delta_ek[0] = -0.35;
    delta_ek[1] = -0.05;
    delta_ek[2] = 0;
    delta_ek[3] = 0.05;
    delta_ek[4] = 0.35;
    fuzzy_control(uk, ek, delta_ek, 2, 100);

    printf("\n");

    return 0;
}