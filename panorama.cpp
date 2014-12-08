/*
 * panorama.cpp
 *
 *  Created on: Mar 28, 2013
 *      Author: Abhiram R (CS10B060)
 *
 */

#include <cv.h>
#include <highgui.h>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <stdlib.h>

using namespace std;

#define MAX_ITER 100
/** Template Code */
#define TINY 1.0e-20;

void lu_decomposition ( double a [9] [9], int n, int indx [], double *d )
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double vv [n + 1];
    *d = 1.0;
    for ( i = 1; i <= n; i++ )
    {
        big = 0.0;
        for ( j = 1; j <= n; j++ )
            if ( (temp = fabs ( a [i] [j] )) > big )
                big = temp;
        if ( big == 0.0 )
            printf ( "Singular matrix in routine ludcmp" );
        vv [i] = 1.0 / big;
    }
    for ( j = 1; j <= n; j++ )
    {
        for ( i = 1; i < j; i++ )
        {
            sum = a [i] [j];
            for ( k = 1; k < i; k++ )
                sum -= a [i] [k] * a [k] [j];
            a [i] [j] = sum;
        }
        big = 0.0;
        for ( i = j; i <= n; i++ )
        {
            sum = a [i] [j];
            for ( k = 1; k < j; k++ )
                sum -= a [i] [k] * a [k] [j];
            a [i] [j] = sum;
            if ( (dum = vv [i] * fabs ( sum )) >= big )
            {
                big = dum;
                imax = i;
            }
        }
        if ( j != imax )
        {
            for ( k = 1; k <= n; k++ )
            {
                dum = a [imax] [k];
                a [imax] [k] = a [j] [k];
                a [j] [k] = dum;
            }
            *d = -(*d);
            vv [imax] = vv [j];
        }
        indx [j] = imax;
        if ( a [j] [j] == 0.0 )
            a [j] [j] = TINY;
        if ( j != n )
        {
            dum = 1.0 / (a [j] [j]);
            for ( i = j + 1; i <= n; i++ )
                a [i] [j] *= dum;
        }
    }
}
void lu_bksb ( double a [9] [9], int n, int indx [], double b [] )
{
    int i, ii = 0, ip, j;
    double sum;

    for ( i = 1; i <= n; i++ )
    {
        ip = indx [i];
        sum = b [ip];
        b [ip] = b [i];
        if ( ii )
            for ( j = ii; j <= i - 1; j++ )
                sum -= a [i] [j] * b [j];
        else if ( sum )
            ii = i;
        b [i] = sum;
    }
    for ( i = n; i >= 1; i-- )
    {
        sum = b [i];
        for ( j = i + 1; j <= n; j++ )
            sum -= a [i] [j] * b [j];
        b [i] = sum / a [i] [i];
    }
}
void ludcmp3 ( double a [4] [4], int n, int indx [], double *d )
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double vv [n + 1];
    *d = 1.0;
    for ( i = 1; i <= n; i++ )
    {
        big = 0.0;
        for ( j = 1; j <= n; j++ )
            if ( (temp = fabs ( a [i] [j] )) > big )
                big = temp;
        vv [i] = 1.0 / big;
    }
    for ( j = 1; j <= n; j++ )
    {
        for ( i = 1; i < j; i++ )
        {
            sum = a [i] [j];
            for ( k = 1; k < i; k++ )
                sum -= a [i] [k] * a [k] [j];
            a [i] [j] = sum;
        }
        big = 0.0;
        for ( i = j; i <= n; i++ )
        {
            sum = a [i] [j];
            for ( k = 1; k < j; k++ )
                sum -= a [i] [k] * a [k] [j];
            a [i] [j] = sum;
            if ( (dum = vv [i] * fabs ( sum )) >= big )
            {
                big = dum;
                imax = i;
            }
        }
        if ( j != imax )
        {
            for ( k = 1; k <= n; k++ )
            {
                dum = a [imax] [k];
                a [imax] [k] = a [j] [k];
                a [j] [k] = dum;
            }
            *d = -(*d);
            vv [imax] = vv [j];
        }
        indx [j] = imax;
        if ( a [j] [j] == 0.0 )
            a [j] [j] = TINY;
        if ( j != n )
        {
            dum = 1.0 / (a [j] [j]);
            for ( i = j + 1; i <= n; i++ )
                a [i] [j] *= dum;
        }
    }
}

void lubksb3 ( double a [4] [4], int n, int indx [], double b [] )
{
    int i, ii = 0, ip, j;
    double sum;

    for ( i = 1; i <= n; i++ )
    {
        ip = indx [i];
        sum = b [ip];
        b [ip] = b [i];
        if ( ii )
            for ( j = ii; j <= i - 1; j++ )
                sum -= a [i] [j] * b [j];
        else if ( sum )
            ii = i;
        b [i] = sum;
    }
    for ( i = n; i >= 1; i-- )
    {
        sum = b [i];
        for ( j = i + 1; j <= n; j++ )
            sum -= a [i] [j] * b [j];
        b [i] = sum / a [i] [i];
    }
}
void ludcmp6 ( float a [7] [7], int n, int indx [], float *d )
{
    int i, imax, j, k;
    float big, dum, sum, temp;
    double vv [n + 1];
    *d = 1.0;
    for ( i = 1; i <= n; i++ )
    {
        big = 0.0;
        for ( j = 1; j <= n; j++ )
            if ( (temp = fabs ( a [i] [j] )) > big )
                big = temp;
        //if ( big == 0.0 )
        //    nrerror ( "Singular matrix in routine ludcmp" );
        vv [i] = 1.0 / big;
    }
    for ( j = 1; j <= n; j++ )
    {
        for ( i = 1; i < j; i++ )
        {
            sum = a [i] [j];
            for ( k = 1; k < i; k++ )
                sum -= a [i] [k] * a [k] [j];
            a [i] [j] = sum;
        }
        big = 0.0;
        for ( i = j; i <= n; i++ )
        {
            sum = a [i] [j];
            for ( k = 1; k < j; k++ )
                sum -= a [i] [k] * a [k] [j];
            a [i] [j] = sum;
            if ( (dum = vv [i] * fabs ( sum )) >= big )
            {
                big = dum;
                imax = i;
            }
        }
        if ( j != imax )
        {
            for ( k = 1; k <= n; k++ )
            {
                dum = a [imax] [k];
                a [imax] [k] = a [j] [k];
                a [j] [k] = dum;
            }
            *d = -(*d);
            vv [imax] = vv [j];
        }
        indx [j] = imax;
        if ( a [j] [j] == 0.0 )
            a [j] [j] = TINY;
        if ( j != n )
        {
            dum = 1.0 / (a [j] [j]);
            for ( i = j + 1; i <= n; i++ )
                a [i] [j] *= dum;
        }
    }
    delete vv;
}

void lubksb6 ( float a [7] [7], int n, int indx [], float b [] )
{
    int i, ii = 0, ip, j;
    float sum;

    for ( i = 1; i <= n; i++ )
    {
        ip = indx [i];
        sum = b [ip];
        b [ip] = b [i];
        if ( ii )
            for ( j = ii; j <= i - 1; j++ )
                sum -= a [i] [j] * b [j];
        else if ( sum )
            ii = i;
        b [i] = sum;
    }
    for ( i = n; i >= 1; i-- )
    {
        sum = b [i];
        for ( j = i + 1; j <= n; j++ )
            sum -= a [i] [j] * b [j];
        b [i] = sum / a [i] [i];
    }
}
void matinv3 ( double** Bin, double** Bout )
{
    double Btemp [4] [4];
    int i, j;
    int indx [4];
    double d, col [4];

    for ( i = 0; i < 3; i++ )
        for ( j = 0; j < 3; j++ )
            Btemp [i + 1] [j + 1] = Bin [i] [j];

    ludcmp3 ( Btemp, 3, indx, &d );
    for ( j = 1; j <= 3; j++ )
    {
        for ( i = 1; i <= 3; i++ )
            col [i] = 0.0;
        col [j] = 1.0;
        lubksb3 ( Btemp, 3, indx, col );
        for ( i = 1; i <= 3; i++ )
            Bout [i - 1] [j - 1] = col [i];
    }
}
/** End of Template Code */
double affine_warp_value ( double** A, int x, int y, int choose )
{
    switch ( choose )
    {
    case 0:
        return (A [0] [0] * x + A [0] [1] * y + A [0] [2]);
    case 1:
        return (A [1] [0] * x + A [1] [1] * y + A [1] [2]);
    default:
        printf ( "error in Affine_warp\n" );
        return (1.0);
    }
}

double warp_value ( double** A, int x, int y, int choose )
{
    switch ( choose )
    {
    case 0:
        return (A [0] [0] * x + A [0] [1] * y + A [0] [2]);
    case 1:
        return (A [1] [0] * x + A [1] [1] * y + A [1] [2]);
    case 2:
        return (A [2] [0] * x + A [2] [1] * y + A [2] [2]);
    default:
        printf ( "error in proj_warp\n" );
        return (1.0);
    }
}

double matching_error ( double** A, int x, int y, IplImage* img_current, IplImage* img_bgd )
{
    //-- Sizes of the images
    CvSize bgd_size = cvGetSize ( img_bgd );
    CvSize cur_size = cvGetSize ( img_current );

    //-- Obtaining the projective warp values
    double denom = warp_value ( A, x, y, 2 );
    double lf_x_prime = warp_value ( A, x, y, 0 ) / denom;
    double lf_y_prime = warp_value ( A, x, y, 1 ) / denom;

    //-- Estimating warp to the nearest pixel
    int x_prime = (int) (lf_x_prime + 0.5);
    int y_prime = (int) (lf_y_prime + 0.5);

    //-- Obtaining the fractional errors
    int frac_x = abs ( lf_x_prime - x_prime );
    int frac_y = abs ( lf_y_prime - y_prime );

    //-- If out of scope, return maximum error: some high value
    if ( x_prime < 0 || y_prime < 0 || x_prime > bgd_size.width || y_prime > bgd_size.height )
        return -1;

    //-- Cyclic restriction
    x_prime = (x_prime + 20 * bgd_size.width) % bgd_size.width;
    y_prime = (y_prime + 20 * bgd_size.height) % bgd_size.height;
    //cout << "xprime = " << x_prime << "; yprime = " << y_prime << endl;

    //-- Error between intensities at x', y' of current image,
    //-- and x, y of the original image
    double e = (1 - frac_x) * (1 - frac_y) * cvGetReal2D ( img_bgd, y_prime, x_prime )
            + frac_x * (1 - frac_y) * cvGetReal2D ( img_bgd, y_prime, (x_prime + 1) % bgd_size.width )
            + (1 - frac_x) * frac_y * cvGetReal2D ( img_bgd, (y_prime + 1) % bgd_size.height, x_prime )
            + frac_x * frac_y * cvGetReal2D ( img_bgd, (y_prime + 1) % bgd_size.height, (x_prime + 1) % bgd_size.width )
            - cvGetReal2D ( img_current, y, x );

    return e;
}
void print_A ( double ** A, int n )
{
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < n; j++ )
            cout << A [i] [j] << " ";
        cout << endl;
    }

}
double sum_squares_error ( IplImage* img_bgd, IplImage* img_current, int x, int y )
{
    double sse = 0;
    for ( int i = 0; i < img_current->width; i++ )
        for ( int j = 0; j < img_current->height; j++ )
        {
            if ( j + y >= img_bgd->height || i + x >= img_bgd->width )
                break;
            double e = cvGetReal2D ( img_current, j, i ) - cvGetReal2D ( img_bgd, j + y, i + x );
            sse += e * e;
        }
    return sse;
}
double** translation_transform ( IplImage* img_bgd, IplImage* img_current, double** A )
{
    CvSize bgd_size = cvGetSize ( img_bgd );
    CvSize cur_size = cvGetSize ( img_current );

    int minx = 0;
    int miny = 0;
    double min_sse = INFINITY;
#pragma omp parallel for schedule(dynamic,1) collapse(2)
    for ( int x = 0; x < img_bgd->width - img_current->width; x++ )
        for ( int y = 0; y < img_bgd->height - img_current->height; y++ )
        {
            double sse = sum_squares_error ( img_bgd, img_current, x, y );
            //cout << sse << endl;
            if ( sse < min_sse )
            {
                minx = x;
                miny = y;
                min_sse = sse;
            }
        }
    A [0] [2] = minx;
    A [1] [2] = miny;
   /* cout << "------------" << endl;
    cout << minx << " " << miny << " " << min_sse << endl;
    cout << "------------" << endl;*/
    return A;

}
double** projective_transform ( IplImage* img_bgd, IplImage* img_current, double** A )
{
    CvSize bgd_size = cvGetSize ( img_bgd );
    CvSize cur_size = cvGetSize ( img_current );

    //-- Normalize A
    for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
            A [i] [j] = A [i] [j] / A [2] [2];

    //-- Optimization matrices for recovering projective transforms
    double optim_lhs [9] [9];
    double optim_rhs [9];
    int optim_index [9];
    double d;

    double pre_sum_errors = 0;
    double sum_errors = 0;

    //-- The Main loop
    int i_loop = 0;
    while ( i_loop < MAX_ITER )
    {
        print_A ( A, 3 );
        pre_sum_errors = 0;
        sum_errors = 0;
        int count_points_1 = 0;
        int count_points_2 = 0;

        //-- Initializing the optim matrices
        for ( int i = 1; i < 9; i++ )
        {
            optim_rhs [i] = 0;
            for ( int j = 1; j < 9; j++ )
                optim_lhs [i] [j] = 0;
        }

        //-- Computing the projective warp of each pixel
//#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for ( int x = 1; x < cur_size.width - 1; x++ )
            for ( int y = 1; y < cur_size.height - 1; y++ )
            {
                double e = matching_error ( A, x, y, img_current, img_bgd );
                if ( e == -1 )
                    continue;
                //cout << " Error in " << x << " " << y << " is " << e << endl;
                pre_sum_errors += e * e;
                count_points_1++;

                //-- Derivates along x and y
                double I_x = (cvGetReal2D ( img_current, y, x + 1 ) - cvGetReal2D ( img_current, y, x - 1 )) / 2;
                double I_y = (cvGetReal2D ( img_current, y + 1, x ) - cvGetReal2D ( img_current, y - 1, x )) / 2;

                //-- Values for updating the warp
                //XXX : Unnecessary re-computation due to function abstraction
                double denom = warp_value ( A, x, y, 2 );
                double lf_x_prime = warp_value ( A, x, y, 0 ) / denom;
                double lf_y_prime = warp_value ( A, x, y, 1 ) / denom;

                double dIdB [9];
                dIdB [1] = I_x * x;
                dIdB [2] = I_x * y;
                dIdB [3] = I_x;
                dIdB [4] = I_y * x;
                dIdB [5] = I_y * y;
                dIdB [6] = I_y;
                dIdB [7] = -x * lf_x_prime * I_x - x * lf_y_prime * I_y;
                dIdB [8] = -y * lf_x_prime * I_x - y * lf_y_prime * I_y;

                //-- Compute A'A and A'e: KLT Estimation of new warp
                for ( int i = 1; i < 9; i++ )
                {
                    optim_rhs [i] -= dIdB [i] * e;
                    for ( int j = 1; j < 9; j++ )
                        optim_lhs [i] [j] += dIdB [i] * dIdB [j];
                }
            }
        //-- End of for loop
        if ( pre_sum_errors != 0 )
        {
            //-- Doing the LU Decomposition
            lu_decomposition ( optim_lhs, 8, optim_index, &d );
            lu_bksb ( optim_lhs, 8, optim_index, optim_rhs );

            //-- Creating the new warp
            double** A_prime = new double* [3];
            for ( int i = 0; i < 3; i++ )
                A_prime [i] = new double [3];
            for ( int i = 1; i < 9; i++ )
            {
                A_prime [(i - 1) / 3] [(i - 1) % 3] = A [(i - 1) / 3] [(i - 1) % 3] + optim_rhs [i];
            }
            A_prime [2] [2] = 1; //XXX : Something is fishy here

            //-- Re-apply the warp and check if it reduces the error

            for ( int x = 1; x < cur_size.width - 1; x++ )
                for ( int y = 1; y < cur_size.height - 1; y++ )
                {
                    double e = matching_error ( A_prime, x, y, img_current, img_bgd );
                    if ( e == -1 )
                        continue;
                    sum_errors += e * e;
                    count_points_2++;

                }
            //-- Update the warp if the error has reduced
            if ( sum_errors / count_points_2 <= pre_sum_errors / count_points_1 )
            {
                for ( int i = 0; i < 3; i++ )
                    for ( int j = 0; j < 3; j++ )
                        A [i] [j] = A_prime [i] [j];
            }
            delete A_prime;
            cout << sum_errors / count_points_2 << " " << pre_sum_errors / count_points_1 << endl;
            //-- If the error actually increases, the exit the loop
            if ( sum_errors / count_points_2 + 0.02 >= pre_sum_errors / count_points_1 )
                i_loop = MAX_ITER;

            i_loop++;
        } else
            i_loop = MAX_ITER;

        //TODO-- Sanity checking XXX ?
    }

}

char* img_1 = "Canon A1200 slow pan test..mp40000.jpg";
char* img_2 = "Canon A1200 slow pan test..mp40001.jpg";
//char* img_1 = "trunks_1.png";
//char* img_2 = "trunks_2.png";
int IMAGES = 50;
char* window_name = "Image Stitching - Abhiram R";
int main ( int argc, char** argv )
{
    IplImage* images [IMAGES];
    for ( int i = 0; i < IMAGES; i++ )
    {

        char* buffer = new char [100];
        sprintf ( buffer, "Canon A1200 slow pan test..mp400%d.jpg", i + 30 );

        images [i] = cvLoadImage ( buffer, CV_LOAD_IMAGE_GRAYSCALE );
    }

    IplImage* trunks_1 = cvCloneImage ( images [0] );

    cvNamedWindow ( window_name, 1 );
    for ( int image_count = 0; image_count < 8; image_count++ )
    {
        IplImage* trunks_2 = cvCloneImage ( images [image_count + 1] );

        cout << cvGetSize ( trunks_1 ).width << endl;
        cout << cvGetSize ( trunks_1 ).height << endl;
        cout << cvGetSize ( trunks_2 ).width << endl;
        cout << cvGetSize ( trunks_2 ).height << endl;

        double** A = new double* [3];
        for ( int i = 0; i < 3; i++ )
            A [i] = new double [3];

        for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 3; j++ )
                if ( i == j )
                    A [i] [j] = 1;
                else
                    A [i] [j] = 0;

        int blur_max = 66;
        bool ssq_estimated = false;
        for ( int i = 5; i < blur_max; i += 10 )
        {
            IplImage* blurred_image_1 = cvCloneImage ( trunks_1 );
            IplImage* blurred_image_2 = cvCloneImage ( trunks_2 );

            cvSmooth ( trunks_1, blurred_image_1, CV_GAUSSIAN, blur_max - i, blur_max - i );
            cvSmooth ( trunks_2, blurred_image_2, CV_GAUSSIAN, blur_max - i, blur_max - i );

            if ( !ssq_estimated )
            {
                ssq_estimated = true;
                translation_transform ( blurred_image_1, blurred_image_2, A );
            }
            projective_transform ( blurred_image_1, blurred_image_2, A );

            cvReleaseImage ( &blurred_image_1 );
            cvReleaseImage ( &blurred_image_2 );
        }
        double** A_p = new double* [3];
        for ( int i = 0; i < 3; i++ )
            A_p [i] = new double [3];
        matinv3 ( A, A_p );
        for ( int i = 0; i < 3; i++ )
        {
            for ( int j = 0; j < 3; j++ )
                cout << A_p [i] [j] << " ";
            cout << endl;
        }

        double minx = INT_MAX;
        double maxx = INT_MIN;
        double miny = INT_MAX;
        double maxy = INT_MIN;
        cout << maxy << endl;
        for ( int x = 1; x <= trunks_2->width; x++ )
            for ( int y = 1; y <= trunks_2->height; y++ )
            {
                double denom = warp_value ( A, x, y, 2 );
                double xx = warp_value ( A, x, y, 0 ) / denom;
                double yy = warp_value ( A, x, y, 1 ) / denom;
                if ( xx < minx )
                    minx = xx;
                if ( xx > maxx )
                    maxx = xx;
                if ( yy < miny )
                    miny = yy;
                if ( yy > maxy )
                    maxy = yy;
            }

        cout << " Params = ";
        cout << minx << " " << maxx << " " << miny << " " << maxy << endl;
        int left_extra = abs ( floor ( min ( minx, 0.0 ) ) );
        int right_extra = ceil ( max ( 0.0, maxx - trunks_1->width ) );
        int top_extra = abs ( floor ( min ( miny, 0.0 ) ) );
        int down_extra = ceil ( max ( 0.0, maxy - trunks_1->height ) );

        cout << left_extra << " " << right_extra << " " << top_extra << " " << down_extra << endl;

        CvSize stitch_size;
        stitch_size.width = max ( trunks_2->width, left_extra + right_extra + trunks_1->width );
        stitch_size.height = max ( trunks_2->height, top_extra + down_extra + trunks_1->height );
        IplImage* stitch_image = cvCreateImage ( stitch_size, trunks_1->depth, trunks_1->nChannels );
        cout << "StitchImage " << stitch_image->width << " " << stitch_image->height << endl;
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for ( int x = 1; x < trunks_1->width; x++ )
            for ( int y = 1; y < trunks_1->height; y++ )
            {
                //cout << x << " " << y << endl;
                cvSet2D ( stitch_image, y + top_extra, x + left_extra, cvGet2D ( trunks_1, y, x ) );

            }
#pragma omp parallel for schedule(dynamic,1) collapse(2)
        for ( int x = 1; x < trunks_2->width; x++ )
            for ( int y = 1; y < trunks_2->height; y++ )
            {
                //cout << x << " " << y <<endl;
                double denom = warp_value ( A, x, y, 2 );
                int xx = (int) (warp_value ( A, x, y, 0 ) / denom) + left_extra; // + 0.5);
                int yy = (int) (warp_value ( A, x, y, 1 ) / denom) + top_extra; // + 0.5);
                //cout << xx << " " << yy <<endl;
                /*if ( ( xx <= left_extra || xx >= stitch_image->width - right_extra )
                 || ( yy <= top_extra
                 || yy >= stitch_image->height - down_extra ) )*/

                cvSet2D ( stitch_image, yy, xx, cvGet2D ( trunks_2, y, x ) );

            }

        cvReleaseImage ( &trunks_1 );
        cvReleaseImage ( &trunks_2 );

        trunks_1 = cvCloneImage ( stitch_image );
        cvReleaseImage ( &stitch_image );

    }
    for ( int i = 0; i < IMAGES; i++ )
    {
        cvReleaseImage ( &images [i] );
    }
    while ( 1 )
    {
        cvShowImage ( window_name, trunks_1 );
        char c = cvWaitKey ( 1 );
        if ( c == 'q' )
            break;
    }

    cvReleaseImage ( &trunks_1 );

}

