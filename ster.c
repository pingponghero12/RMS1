// ster.c Manfred Gawlas

#include <stdio.h>
#include <math.h>

struct AbGeometry // Aktualizowane rozmiary ziarna
{
	double L;
	double R;
	double D;
};

struct General // Zbiór wielkości które przekazujemy do funkcji a nie wsadziłęm ich gdzie indziej
{
	double Pch;
	double L;
	double D;
	double d;
	double At;
	double Pe;
};

struct Fuel // Wielkości paliwa
{
	double Cstar;
	double a;
	double n;
	double k;
	double rho;
};

double CF(struct General Cylinder, struct Fuel RNX71) // Funkcja zwracająca wartość CF
{
	return pow((2*RNX71.k*RNX71.k)/(RNX71.k-1) * pow(2/(RNX71.k+1), (RNX71.k+1)/(RNX71.k-1)) * (1 - pow(Cylinder.Pe/Cylinder.Pch, (RNX71.k-1)/RNX71.k)), 0.5);
}

double F(struct General Cylinder, struct Fuel RNX71) // Fukcja zwracająca wartość ciągu
{
	return CF(Cylinder, RNX71) * Cylinder.At * Cylinder.Pch;
}

double FunctionAb(struct AbGeometry Ab) // Funkcja zwracająca area burning
{
	return 3.14159*(2*Ab.L*Ab.R + 3*0.25*Ab.D*Ab.D - 3*Ab.R*Ab.R);
}


int main(void)
{
	// Wszystkie wartości są wyrażone w podstawowych jednostkach układu SI, to jest ciśnienie [Pa], długości [m], pola [m^2], itd.

	double r=0;
	double Ic=0;
	double F0=82.37528;
	double tc=0;
	double totalRegressed=0; // Zmienna przechowywujące ile do tej pory długości w głąb spaliło się [m]
	double Deltat=0.003; // Zmiana czasu dla ktorej liczymy prostokąty

	struct General Cylinder;
	Cylinder.Pch=35*101325; 
	Cylinder.L=0.15;
	Cylinder.D=0.023;
	Cylinder.d=0.01;
	//Cylinder.At=0.00001824875;
	Cylinder.Pe=101325;
	
	struct Fuel RNX71;
	RNX71.Cstar=779;
	RNX71.a=0.0000163938; // Współczynnik dla Pa, wyliczony z wsp. dla MPa
	RNX71.n=0.371;
	RNX71.k=1.18;
	RNX71.rho=1848;

	struct AbGeometry Ab;
	Ab.L=Cylinder.L;
	Ab.R=Cylinder.d/2;
	Ab.D=Cylinder.D;
	
	//printf("F0=%f\n", F(Cylinder, RNX71));
	//printf("CF=%f\n\n", CF(Cylinder, RNX71));
	

	Cylinder.At=F0/(Cylinder.Pch*CF(Cylinder, RNX71));
	Cylinder.Pch=pow(FunctionAb(Ab) / Cylinder.At, 1/(1-RNX71.n)) * pow(RNX71.Cstar * RNX71.rho * RNX71.a, 1/(1-RNX71.n));
	
	//printf("de = %f\n", 2 * sqrt(Ae / 3.1415));
	printf("dt = %f\n", 2 * sqrt(Cylinder.At / 3.1415));
	printf("F0=%f\n", F(Cylinder, RNX71));
	//////////////////////////////////
	// Główna pętla programu, algorytm podobny do metody liczenia całki przez prostokąty. Opisany szczegółowo w dokumentacji
	//////////////////////////////////
	
	while(totalRegressed<((Cylinder.D-Cylinder.d)/2))
	{
		// Obliczanie wartości ciągu oraz Impulsu, liczenie czasu

		tc+=Deltat;
		F0=F(Cylinder, RNX71);
		Ic+=F0 * Deltat;
		
		
		///////////////////////////////
		// Poniżej miejsce na printowanie wartości dla każdej iteracji pętli. Pozwala zbierać dane do wykresów.
		///////////////////////////////
		
		//printf("%f\n", F0);
		//printf("%f\n\n", Ic);
		//printf("%f\n", FunctionAb(Ab)/Cylinder.At); // Kn
	//	printf("%f\n", Cylinder.Pch);
		
		///////////////////////////////

		r=RNX71.a * pow(Cylinder.Pch, RNX71.n);
		Ab.L=Ab.L - 3*r*Deltat;
		Ab.R=Ab.R + r*Deltat;

		Cylinder.Pch=pow(FunctionAb(Ab) / Cylinder.At, 1/(1-RNX71.n)) * pow(RNX71.Cstar * RNX71.rho * RNX71.a, 1/(1-RNX71.n));
		
		totalRegressed+=r*Deltat; // Uaktualnianie wartości spalonej części do tej pory
	}


	// Printowanie wartości końcowych oraz impulsu
	
	printf("Ic=%f\n", Ic);
	printf("tc=%f\n", tc);	
	printf("Pch=%f\n", Cylinder.Pch);	
	printf("At=%f\n", Cylinder.At);

	return 0;
}
