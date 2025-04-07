    int Nmax = 100000;
    gridGenerate();
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < N - 1; j++) {
            cout << Gij[i][j] << "   ";
        }
        cout << endl;
    }
    double a = solveLaplaceEquationAtOnePoint(2, 2, Nmax);
    cout << "Result: " << a << endl;
    return 0;