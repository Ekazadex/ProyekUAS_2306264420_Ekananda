#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <functional>
#include <string>
#include <algorithm>

using namespace std;
using namespace std::chrono;

class RungeKuttaSolver {
private:
    double h;           // step size
    double x0, y0;      // initial conditions
    double xn;          // end point
    function<double(double, double)> f;  // ODE function
    function<double(double)> analytical; // analytical solution (for comparison)
    
public:
    RungeKuttaSolver(double step, double x_init, double y_init, double x_end, 
                     function<double(double, double)> func,
                     function<double(double)> analytical_func = nullptr) 
        : h(step), x0(x_init), y0(y_init), xn(x_end), f(func), analytical(analytical_func) {}
    
    vector<pair<double, double>> solve();
    double calculateMaxError();
    double calculateRelativeError();
    void printResults();
    void saveToFile(const string& filename);
    double getComputationTime();
    void performConvergenceAnalysis();
    void compareWithOtherMethods();
    
    // Additional methods for other numerical methods
    vector<pair<double, double>> solveEuler();
    vector<pair<double, double>> solveRK2();
    vector<pair<double, double>> solveRK3();
};

vector<pair<double, double>> RungeKuttaSolver::solve() {
    vector<pair<double, double>> result;
    double x = x0, y = y0;
    
    result.push_back(make_pair(x, y));
    
    while (x < xn - h/2) {  // Avoid floating point precision issues
        // Runge-Kutta 4th order method implementation
        double k1 = h * f(x, y);
        double k2 = h * f(x + h/2.0, y + k1/2.0);
        double k3 = h * f(x + h/2.0, y + k2/2.0);
        double k4 = h * f(x + h, y + k3);
        
        // RK4 formula
        y = y + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
        x = x + h;
        
        result.push_back(make_pair(x, y));
    }
    
    return result;
}

vector<pair<double, double>> RungeKuttaSolver::solveEuler() {
    vector<pair<double, double>> result;
    double x = x0, y = y0;
    
    result.push_back(make_pair(x, y));
    
    while (x < xn - h/2) {
        // Euler's method: y_{n+1} = y_n + h*f(x_n, y_n)
        y = y + h * f(x, y);
        x = x + h;
        result.push_back(make_pair(x, y));
    }
    
    return result;
}

vector<pair<double, double>> RungeKuttaSolver::solveRK2() {
    vector<pair<double, double>> result;
    double x = x0, y = y0;
    
    result.push_back(make_pair(x, y));
    
    while (x < xn - h/2) {
        // RK2 (Heun's method)
        double k1 = h * f(x, y);
        double k2 = h * f(x + h, y + k1);
        
        y = y + (k1 + k2) / 2.0;
        x = x + h;
        result.push_back(make_pair(x, y));
    }
    
    return result;
}

vector<pair<double, double>> RungeKuttaSolver::solveRK3() {
    vector<pair<double, double>> result;
    double x = x0, y = y0;
    
    result.push_back(make_pair(x, y));
    
    while (x < xn - h/2) {
        // RK3 method
        double k1 = h * f(x, y);
        double k2 = h * f(x + h/2.0, y + k1/2.0);
        double k3 = h * f(x + h, y - k1 + 2.0*k2);
        
        y = y + (k1 + 4.0*k2 + k3) / 6.0;
        x = x + h;
        result.push_back(make_pair(x, y));
    }
    
    return result;
}

double RungeKuttaSolver::calculateMaxError() {
    if (!analytical) return -1.0;  // No analytical solution provided
    
    auto results = solve();
    double maxError = 0.0;
    
    for (const auto& point : results) {
        double x = point.first;
        double y_numerical = point.second;
        double y_analytical = analytical(x);
        double error = abs(y_numerical - y_analytical);
        if (error > maxError) {
            maxError = error;
        }
    }
    
    return maxError;
}

double RungeKuttaSolver::calculateRelativeError() {
    if (!analytical) return -1.0;
    
    auto results = solve();
    double maxRelError = 0.0;
    
    for (const auto& point : results) {
        double x = point.first;
        double y_numerical = point.second;
        double y_analytical = analytical(x);
        if (abs(y_analytical) > 1e-10) {  // Avoid division by very small numbers
            double relError = abs((y_numerical - y_analytical) / y_analytical);
            if (relError > maxRelError) {
                maxRelError = relError;
            }
        }
    }
    
    return maxRelError;
}

void RungeKuttaSolver::printResults() {
    auto results = solve();
    
    cout << fixed << setprecision(8);
    cout << "\n=== HASIL METODE RUNGE-KUTTA ORDE 4 ===\n";
    cout << "Step size (h): " << h << "\n";
    cout << "Interval: [" << x0 << ", " << xn << "]\n";
    cout << "Kondisi awal: y(" << x0 << ") = " << y0 << "\n";
    cout << "Jumlah langkah: " << results.size() - 1 << "\n\n";
    
    cout << setw(12) << "x" << setw(16) << "y (RK4)";
    if (analytical) {
        cout << setw(16) << "y (Analytical)" << setw(16) << "Absolute Error" << setw(16) << "Relative Error";
    }
    cout << "\n" << string(analytical ? 76 : 28, '-') << "\n";
    
    for (const auto& point : results) {
        double x = point.first;
        double y_num = point.second;
        
        cout << setw(12) << x << setw(16) << y_num;
        
        if (analytical) {
            double y_anal = analytical(x);
            double absError = abs(y_num - y_anal);
            double relError = (abs(y_anal) > 1e-10) ? abs((y_num - y_anal) / y_anal) : 0.0;
            cout << setw(16) << y_anal << setw(16) << scientific << absError 
                 << setw(16) << relError << fixed;
        }
        cout << "\n";
    }
    
    if (analytical) {
        cout << "\n=== ANALISIS ERROR ===\n";
        cout << "Max Absolute Error: " << scientific << calculateMaxError() << "\n";
        cout << "Max Relative Error: " << calculateRelativeError() << "\n";
        cout << fixed;
    }
}

void RungeKuttaSolver::saveToFile(const string& filename) {
    auto results = solve();
    ofstream file(filename);
    
    file << "x,y_rk4";
    if (analytical) file << ",y_analytical,absolute_error,relative_error";
    file << "\n";
    
    for (const auto& point : results) {
        double x = point.first;
        double y_num = point.second;
        
        file << fixed << setprecision(10) << x << "," << y_num;
        
        if (analytical) {
            double y_anal = analytical(x);
            double absError = abs(y_num - y_anal);
            double relError = (abs(y_anal) > 1e-10) ? abs((y_num - y_anal) / y_anal) : 0.0;
            file << "," << y_anal << "," << scientific << absError << "," << relError << fixed;
        }
        file << "\n";
    }
    
    file.close();
}

double RungeKuttaSolver::getComputationTime() {
    auto start = high_resolution_clock::now();
    solve();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    return duration.count() / 1000000.0;  // Convert to seconds
}

void RungeKuttaSolver::performConvergenceAnalysis() {
    if (!analytical) {
        cout << "Error: Analytical solution required for convergence analysis!\n";
        return;
    }
    
    cout << "\n=== ANALISIS KONVERGENSI ===\n";
    cout << setw(15) << "Step Size (h)" << setw(15) << "Error Max" << setw(10) << "Order" << setw(15) << "CPU Time (s)" << "\n";
    cout << string(55, '-') << "\n";
    
    vector<double> step_sizes = {0.2, 0.1, 0.05, 0.025};
    vector<double> errors;
    
    for (double step : step_sizes) {
        RungeKuttaSolver temp_solver(step, x0, y0, xn, f, analytical);
        double error = temp_solver.calculateMaxError();
        double comp_time = temp_solver.getComputationTime();
        errors.push_back(error);
        
        double order = -1.0;
        if (errors.size() > 1) {
            order = log(errors[errors.size()-2] / error) / log(2.0);
        }
        
        cout << setw(15) << step << setw(15) << scientific << error;
        if (order > 0) {
            cout << setw(10) << fixed << setprecision(1) << order;
        } else {
            cout << setw(10) << "-";
        }
        cout << setw(15) << fixed << setprecision(6) << comp_time << "\n";
    }
}

void RungeKuttaSolver::compareWithOtherMethods() {
    if (!analytical) {
        cout << "Error: Analytical solution required for method comparison!\n";
        return;
    }
    
    cout << "\n=== PERBANDINGAN METODE NUMERIK ===\n";
    cout << setw(15) << "Metode" << setw(10) << "Order" << setw(15) << "Error Max" << setw(20) << "Evaluasi f per step" << "\n";
    cout << string(60, '-') << "\n";
    
    // Euler method
    auto euler_results = solveEuler();
    double euler_error = 0.0;
    for (const auto& point : euler_results) {
        double x = point.first;
        double y_num = point.second;
        double y_anal = analytical(x);
        double error = abs(y_num - y_anal);
        if (error > euler_error) euler_error = error;
    }
    
    // RK2 method
    auto rk2_results = solveRK2();
    double rk2_error = 0.0;
    for (const auto& point : rk2_results) {
        double x = point.first;
        double y_num = point.second;
        double y_anal = analytical(x);
        double error = abs(y_num - y_anal);
        if (error > rk2_error) rk2_error = error;
    }
    
    // RK3 method
    auto rk3_results = solveRK3();
    double rk3_error = 0.0;
    for (const auto& point : rk3_results) {
        double x = point.first;
        double y_num = point.second;
        double y_anal = analytical(x);
        double error = abs(y_num - y_anal);
        if (error > rk3_error) rk3_error = error;
    }
    
    // RK4 method
    double rk4_error = calculateMaxError();
    
    cout << setw(15) << "Euler" << setw(10) << "1" << setw(15) << scientific << euler_error << setw(20) << "1" << "\n";
    cout << setw(15) << "RK2 (Heun)" << setw(10) << "2" << setw(15) << scientific << rk2_error << setw(20) << "2" << "\n";
    cout << setw(15) << "RK3" << setw(10) << "3" << setw(15) << scientific << rk3_error << setw(20) << "3" << "\n";
    cout << setw(15) << "RK4" << setw(10) << "4" << setw(15) << scientific << rk4_error << setw(20) << "4" << "\n";
}

// Test Functions Namespace
namespace TestFunctions {
    // Test Case 1: dy/dx = y, y(0) = 1, analytical: y = e^x
    double case1_ode(double x, double y) {
        return y;
    }
    
    double case1_analytical(double x) {
        return exp(x);
    }
    
    // Test Case 2: dy/dx = -2xy, y(0) = 1, analytical: y = e^(-x²)
    double case2_ode(double x, double y) {
        return -2.0 * x * y;
    }
    
    double case2_analytical(double x) {
        return exp(-x*x);
    }
    
    // Test Case 3: dy/dx = x + y, y(0) = 0, analytical: y = e^x - x - 1
    double case3_ode(double x, double y) {
        return x + y;
    }
    
    double case3_analytical(double x) {
        return exp(x) - x - 1.0;
    }
}

// Analysis Tools
namespace AnalysisTools {
    void testStability(function<double(double, double)> f, double lambda, double h_max) {
        cout << "\n=== UJI STABILITAS ===\n";
        cout << "Persamaan uji: dy/dx = " << lambda << "*y\n";
        cout << "Testing stability dengan berbagai step sizes...\n\n";
        
        vector<double> test_h = {h_max, h_max/2, h_max/4, h_max/8};
        
        for (double h : test_h) {
            RungeKuttaSolver solver(h, 0.0, 1.0, 2.0, 
                [lambda](double x, double y) { return lambda * y; },
                [lambda](double x) { return exp(lambda * x); });
            
            auto results = solver.solve();
            bool stable = true;
            
            for (const auto& point : results) {
                if (abs(point.second) > 1e6) {  // Check for numerical blow-up
                    stable = false;
                    break;
                }
            }
            
            cout << "h = " << h << ": " << (stable ? "STABIL" : "TIDAK STABIL") << "\n";
        }
    }
}

// Output Manager
namespace OutputManager {
    void generateFullReport(RungeKuttaSolver& solver, const string& case_name) {
        cout << "\n" << string(80, '=') << "\n";
        cout << "LAPORAN LENGKAP - " << case_name << "\n";
        cout << string(80, '=') << "\n";
        
        // Print basic results
        solver.printResults();
        
        // Performance analysis
        cout << "\n=== ANALISIS PERFORMA ===\n";
        cout << "Waktu komputasi: " << fixed << setprecision(6) << solver.getComputationTime() << " detik\n";
        
        // Convergence analysis
        solver.performConvergenceAnalysis();
        
        // Method comparison
        solver.compareWithOtherMethods();
        
        // Save to file
        string filename = case_name + "_results.csv";
        solver.saveToFile(filename);
        cout << "\nHasil disimpan ke file: " << filename << "\n";
    }
}

int main() {
    cout << "=== IMPLEMENTASI METODE RUNGE-KUTTA ORDE 4 ===\n";
    cout << "Proyek UAS Metode Numerik\n";
    cout << "Oleh: Ekananda Zhafif Dean\n\n";
    
    // Test Case 1: Exponential Growth
    cout << "\n" << string(80, '=') << "\n";
    cout << "KASUS UJI 1: PERTUMBUHAN EKSPONENSIAL\n";
    cout << "dy/dx = y, y(0) = 1\n";
    cout << "Solusi analitik: y(x) = e^x\n";
    cout << string(80, '=') << "\n";
    
    RungeKuttaSolver case1(0.1, 0.0, 1.0, 1.0, 
                          TestFunctions::case1_ode, 
                          TestFunctions::case1_analytical);
    
    OutputManager::generateFullReport(case1, "Kasus_1_Eksponensial");
    
    // Test Case 2: Gaussian Decay
    cout << "\n" << string(80, '=') << "\n";
    cout << "KASUS UJI 2: PELURUHAN GAUSSIAN\n";
    cout << "dy/dx = -2xy, y(0) = 1\n";
    cout << "Solusi analitik: y(x) = e^(-x²)\n";
    cout << string(80, '=') << "\n";
    
    RungeKuttaSolver case2(0.1, 0.0, 1.0, 1.0, 
                          TestFunctions::case2_ode, 
                          TestFunctions::case2_analytical);
    
    OutputManager::generateFullReport(case2, "Kasus_2_Gaussian");
    
    // Test Case 3: Linear Non-homogeneous
    cout << "\n" << string(80, '=') << "\n";
    cout << "KASUS UJI 3: PERSAMAAN LINEAR NON-HOMOGEN\n";
    cout << "dy/dx = x + y, y(0) = 0\n";
    cout << "Solusi analitik: y(x) = e^x - x - 1\n";
    cout << string(80, '=') << "\n";
    
    RungeKuttaSolver case3(0.1, 0.0, 0.0, 1.0, 
                          TestFunctions::case3_ode, 
                          TestFunctions::case3_analytical);
    
    OutputManager::generateFullReport(case3, "Kasus_3_Linear");
    
    // Stability Analysis
    cout << "\n" << string(80, '=') << "\n";
    cout << "ANALISIS STABILITAS\n";
    cout << string(80, '=') << "\n";
    
    AnalysisTools::testStability(
        [](double x, double y) { return -1.0 * y; }, -1.0, 2.8);
    
    AnalysisTools::testStability(
        [](double x, double y) { return -10.0 * y; }, -10.0, 0.28);
    
    AnalysisTools::testStability(
        [](double x, double y) { return -100.0 * y; }, -100.0, 0.028);
    
    cout << "\n" << string(80, '=') << "\n";
    cout << "ANALISIS SELESAI\n";
    cout << "Semua hasil telah disimpan dalam file CSV\n";
    cout << string(80, '=') << "\n";
    
    return 0;
}