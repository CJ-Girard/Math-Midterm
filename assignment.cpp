#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

// First function: the formula equal to the derivative of y with respect to t
long double f1(long double t, long double y){
    return (t*exp(-y))+(t/(1+pow(t,2)));
}

// Generates the solution to the first IVP
long double y1(long double t){
    return (log(1+pow(t,2)));
}

// Second function: the formula for the secont IVP
long double f2(long double t, long double y){
    long double a = y/2;
    long double b = (3.0/2.0)*(t-1);
    return (a+b);
}

// Generates the solution to the second IVP
long double y2(long double t){
    return (4*exp(t/2)-3*t-3);
}

// Can return either one depending on parameter
long double f(int i, long double t, long double w){
    // Uses the i variable to determine if we're using f1 or f2:
    switch (i) {
    case 1:
        return f1(t, w);
        break;
    case 2:
        return f2(t, w);
        break;
    default:
        cout << "Error in passing function type to f\n"; 
        return -1;
    }
}

// Can return either y depending on parameter
long double findy(int i, long double t){
    switch (i) {
    case 1:
        return y1(t);
        break;
    case 2:
        return y2(t);
        break;
    default:
        cout << "Error in passing function type to y\n"; 
        return -1;
    }
}

// Performs the modified euler method
long double modEuler(int i, long double h, long double t, long double w){
    long double k1=f(i, t, w);
    long double k2=f(i, (t+h), (w+h*k1));
    return w+((h/2)*(k1+k2));
}

void modEulerTest(int i, int N){
    cout << "\n\n\nModified Euler on Formula 1 -- " << N << " iterations:\n";
    cout << "Test:   Formula:      Steps:      Iteration j:           w[j]:               y(t[j]):            Error:" << endl;
    long double w[N+1];
    long double t[N+1];
    long double y[N+1];
    //modEuler will populate t and w matrices
    long double Ndouble = N; //Needs to become a double in order to be passed on into the modEuluer function.
    
    //This part was originally in the modEuler function, but it shouldn be here instead:
    // This should always be used for each of the methods; sets interval, initial values, and h size depending on which case we want to test.
    long double A;
    long double B;
    long double h=1/Ndouble;
    //Switch function sets initial values based on the test we're running:
    switch(i){
    case 1:
        A=0;
        B=2;
        w[0]=0;
        break;
    case 2:
        A=0;
        B=3;
        w[0]=1;
        break;
    default:
        cout << "Error in passing function type to k1\n";
        return;
    }
    h*=(B-A);
    t[0]=A;
    
    for(int j=0; j<N; j++){
        t[j+1]=t[j]+h;
        w[j+1]=modEuler(i, h, t[j], w[j]);
    }

    //Another for loop to populate y values
    for(int j=0; j<=N; j++){
        y[j]=findy(i, t[j]);
    }
    for(int j=0; j<=N; j++){
        cout << "Mod-E    " << i << "   " << N << "   " << j << "   " << w[j] << "   " << y[j] << "   " << abs(y[j]-w[j]) << endl;
    }
    return;
}

long double midpoint(int i, long double h, long double t, long double w){
    long double k1=f(i, t, w);
    long double k2=f(i, (t+(h/2)), (w+(h/2)*k1));
    return w+h*k2;
}

void midpointTest(int i, int N){
    cout << "\n\n\nMidpoit on Formula 1 -- " << N << " iterations:\n";
    cout << "Test:   Formula:      Steps:      Iteration j:           w[j]:               y(t[j]):            Error:" << endl;
    long double w[N+1];
    long double t[N+1];
    long double y[N+1];
    //modEuler will populate t and w matrices
    long double Ndouble = N; //Needs to become a double in order to be passed on into the modEuluer function.
    // This should always be used for each of the methods; sets interval, initial values, and h size depending on which case we want to test.
    long double A;
    long double B;
    long double h=1/Ndouble;
    //Switch function sets initial values based on the test we're running:
    switch(i){
    case 1:
        A=0;
        B=2;
        w[0]=0;
        break;
    case 2:
        A=0;
        B=3;
        w[0]=1;
        break;
    default:
        cout << "Error in passing function type to k1\n";
        return;
    }
    h*=(B-A);
    t[0]=A;
    
    for(int j=0; j<N; j++){
        t[j+1]=t[j]+h;
        w[j+1]=midpoint(i, h, t[j], w[j]);
    }
    //Another for loop to populate y values
    for(int j=0; j<=N; j++){
        y[j]=findy(i, t[j]);
    }
    for(int j=0; j<=N; j++){
        cout << "Midpoint    " << i << "   " << N << "   " << j << "   " << w[j] << "   " << y[j] << "   " << abs(y[j]-w[j]) << endl;
    }
    return;
}

long double fourthOrder(int i, long double h, long double t, long double w){
    //Needs to initialize a for loop that fills every value in the array. We already have w[0] so next we need to fill w[1], w[2], ... , w[N].
        long double k1=f(i, t, w);
        long double k2=f(i, (t+(h/2)), (w+(h/2)*k1));
        long double k3=f(i, t+(h/2), (w+(h/2)*k2));
        long double k4=f(i, t+h, w+h*k3);
        return w+(h/6)*(k1+2*k2+2*k3+k4);
}

void fourthOrderTest(int i, int N){
    cout << "\n\n\nFourth Order on Formula 1 -- " << N << " iterations:\n";
    cout << "Test:   Formula:      Steps:      Iteration j:           w[j]:               y(t[j]):            Error:" << endl;
    long double w[N+1];
    long double t[N+1];
    long double y[N+1];
    //fourthOrder will populate t and w matrices
    long double Ndouble = N; //Needs to become a double in order to be passed on into the modEuluer function.
    // This should always be used for each of the methods; sets interval, initial values, and h size depending on which case we want to test.
    long double A;
    long double B;
    long double h=1/Ndouble;
    //Switch function sets initial values based on the test we're running:
    switch(i){
    case 1:
        A=0;
        B=2;
        w[0]=0;
        break;
    case 2:
        A=0;
        B=3;
        w[0]=1;
        break;
    default:
        cout << "Error in passing function type to k1\n";
        return;
    }
    h*=(B-A);
    t[0]=A;
    for(int j=0; j<N; j++){
        t[j+1]=t[j]+h;
        w[j+1]=fourthOrder(i, h, t[j], w[j]);
    }
    //Another for loop to populate y values
    for(int j=0; j<=N; j++){
        y[j]=findy(i, t[j]);
    }
    for(int j=0; j<=N; j++){
        cout << "Fourth Order:    " << i << "   " << N << "   " << j << "   " << w[j] << "   " << y[j] << "   " << abs(y[j]-w[j]) << endl;
    }
    return;
}

long double richardson(int startOrder, int iterations, long double Nin[]){
    if(iterations==0){return Nin[0];}   //This is a recursive function, and will return itself unless it passes the information to stop, which is when iterations == 0.
    else{
        long double Nout[iterations];   //If it's the last iteration (iterations == 1), only one output as the answer. Second to last, 2 outputs for the remaining iteration. Third to last, 3 outputs, etc ...
        for(int j=0; j<iterations; j++){
            Nout[j]=(pow(2,startOrder)*Nin[j+1]-Nin[j])/(pow(2,startOrder)-1);  //Each Nin[j] is N(h), N(h/2), N(h/4), ...
        }
        startOrder+=1;  //The next iteration will worry about the next power in the sequence.
        iterations-=1;  //The next iteration will consider itself to be the first iteration, only starting from a different power and with one less iteration to go before stopping.
        return richardson(startOrder, iterations, Nout);
    }
}

void modEulerRichLocal(int i, int N, int iterations){
    long double w[N+1];
    long double t[N+1];
    long double y[N+1];
    long double Ndouble = N;
    long double A;
    long double B;
    long double h=1/Ndouble;
    switch(i){
    case 1:
        A=0;
        B=2;
        w[0]=0;
        break;
    case 2:
        A=0;
        B=3;
        w[0]=1;
        break;
    default:
        cout << "Error in passing function type to k1\n";
        return;
    }
    h*=(B-A);
    t[0]=A;
    
    for(int j=0; j<N; j++){
        t[j+1]=t[j]+h;
        long double Nin[iterations+1];  //Array size is iterations+1 because for one iteration, requires 2 N1's, 2 iterations requires 3, ...
        for(int k=0; k<=iterations; k++){   //Will loop to populate all the N1s for h, h/2, h/4, .... where k represents the power of 2 for this particular loop. Note that it will run for the exact no. of times it needs to fill Nin[].
            int smallwindex=pow(2,k)+1;
            long double _h = h/pow(2,k);   //_h is the shortened step size, like 
            long double _w[smallwindex];    //the temporary array of intermediate values for w. smallwindex is the number of steps required to get to t(j+1) and including t[0] being the initial value for t(j) we're left with this index size.
            long double _t[smallwindex];    //the temporary array of intermediate values for t
            _w[0]=w[j];
            _t[0]=t[j];
            for(int a=0; a<pow(2,k); a++){
                _t[a+1]=_t[a]+_h;
                _w[a+1]=modEuler(i, _h, _t[a], _w[a]);
            }
            Nin[k]=_w[smallwindex-1];
        }
        w[j+1]=richardson(3, iterations, Nin);  //Modified Euler and Midpoint have to start at 3 because the first order of h to cancel is O(h^3)
    }
    //Another for loop to populate y values
    for(int j=0; j<=N; j++){
        y[j]=findy(i, t[j]);
    }
    for(int j=0; j<=N; j++){
        cout << "Mod-E Richardson   " << i << "   " << N << "   " << j << "   " << w[j] << "   " << y[j] << "   " << abs(y[j]-w[j]) << endl;
    }
    return;
}

void midpointRichLocal(int i, int N, int iterations){
    long double w[N+1];
    long double t[N+1];
    long double y[N+1];
    long double Ndouble = N;
    long double A;
    long double B;
    long double h=1/Ndouble;
    switch(i){
    case 1:
        A=0;
        B=2;
        w[0]=0;
        break;
    case 2:
        A=0;
        B=3;
        w[0]=1;
        break;
    default:
        cout << "Error in passing function type to k1\n";
        return;
    }
    h*=(B-A);
    t[0]=A;
    
    for(int j=0; j<N; j++){
        t[j+1]=t[j]+h;
        long double Nin[iterations+1];
        for(int k=0; k<=iterations; k++){
            int smallwindex=pow(2,k)+1;
            long double _h = h/pow(2,k);
            long double _w[smallwindex];
            long double _t[smallwindex];
            _w[0]=w[j];
            _t[0]=t[j];
            for(int a=0; a<pow(2,k); a++){
                _t[a+1]=_t[a]+_h;
                _w[a+1]=midpoint(i, _h, _t[a], _w[a]);
            }
            Nin[k]=_w[smallwindex-1];
        }
        w[j+1]=richardson(3, iterations, Nin);
    }
    for(int j=0; j<=N; j++){
        y[j]=findy(i, t[j]);
    }
    for(int j=0; j<=N; j++){
        cout << "Midpoint Richardson   " << i << "   " << N << "   " << j << "   " << w[j] << "   " << y[j] << "   " << abs(y[j]-w[j]) << endl;
    }
    return;
}

void fourthOrderRichLocal(int i, int N, int iterations){
    long double w[N+1];
    long double t[N+1];
    long double y[N+1];
    long double Ndouble = N;
    long double A;
    long double B;
    long double h=1/Ndouble;
    switch(i){
    case 1:
        A=0;
        B=2;
        w[0]=0;
        break;
    case 2:
        A=0;
        B=3;
        w[0]=1;
        break;
    default:
        cout << "Error in passing function type to k1\n";
        return;
    }
    h*=(B-A);
    t[0]=A;
    
    for(int j=0; j<N; j++){
        t[j+1]=t[j]+h;
        long double Nin[iterations+1];
        for(int k=0; k<=iterations; k++){
            int smallwindex=pow(2,k)+1;
            long double _h = h/pow(2,k);
            long double _w[smallwindex];
            long double _t[smallwindex];
            _w[0]=w[j];
            _t[0]=t[j];
            for(int a=0; a<pow(2,k); a++){
                _t[a+1]=_t[a]+_h;
                _w[a+1]=fourthOrder(i, _h, _t[a], _w[a]);
            }
            Nin[k]=_w[smallwindex-1];
        }
        w[j+1]=richardson(5, iterations, Nin);
    }
    for(int j=0; j<=N; j++){
        y[j]=findy(i, t[j]);
    }
    for(int j=0; j<=N; j++){
        cout << "Fourth Order   " << i << "   " << N << "   " << j << "   " << w[j] << "   " << y[j] << "   " << abs(y[j]-w[j]) << endl;
    }
    return;
}

void modEulerRichGlobal(int i, int N, int iterations){
    long double Nin[iterations+1];
    for(int a=0; a<iterations+1; )
}



int main(){
    cout.precision(20);
    fourthOrderRichLocal(2, 64, 12);
    return 0;
}
