// C libraries
#include <stdio.h> // printf()
#include <time.h> // time()
#include <stdlib.h> // srand()

// C++ libraries
#include <string> // string
#include <iostream> // cout, endl
#include <bitset> // bitset
#include <vector> // vector
#include <algorithm> // min(), max()

using std::string;
using std::cout;
using std::endl;
using std::bitset;
using std::vector;
using std::min;
using std::max;

class RSASystem{
public:
	// trace of the last run of millerRabin_verbose()
	vector<int> millerRabinTrace;

	RSASystem(){
		// random number seed
		srand(time(NULL));
	}

	int getn(){
		int answer = 0b1000001;

		// gets the five desired bits
		for (int i = 1 ; i <= 5; i++) {

			// creates the random int
			int randomInt = rand();

			// extracts the right most bit
			int rightMostBit = randomInt & 0b1;

			// puts the desired bit in the correct place
			int desiredPlace = (rightMostBit << 1) & -2 & 2;
			desiredPlace = desiredPlace << (i-1);

			// adds the desired bit to the appropriate place in answer
			answer += desiredPlace;
		}
		return answer;
	}

	int getn_verbose(){
		int answer = 0b1000001;

		// gets the five desired bits
		for (int i = 1 ; i <= 5; i++) {

			// creates the random int
			int randomInt = rand();
			cout << "Random int: " << (bitset<32>) randomInt << endl;

			// extracts the right most bit
			int rightMostBit = randomInt & 0b1;
			cout << "Extracted bit: " << (bitset<1>) rightMostBit << endl;

			// puts the desired bit in the correct place
			int desiredPlace = (rightMostBit << 1) & -2 & 2;
			desiredPlace = desiredPlace << (i-1);

			// adds the desired bit to the appropriate place in answer
			answer += desiredPlace;

			printf("\n");
		}
		return answer;
	}

	int geta(int n){
		int answer = rand()%n;
		while (answer == 0) {
			answer = rand()%n;
		}
		return answer;
	}

	vector<int> getx(int n){
		vector<int> x;
		// gets bits 0 to 6 of n and puts them in x
		int temp = n - 1;
		x.push_back(temp & 0b1);
		x.push_back((temp & 0b10) >> 1);
		x.push_back((temp & 0b100) >> 2);
		x.push_back((temp & 0b1000) >> 3);
		x.push_back((temp & 0b10000) >> 4);
		x.push_back((temp & 0b100000) >> 5);
		x.push_back((temp & 0b1000000) >> 6);
		return x;
	}

	// performs a Miller Rabin test
	bool millerRabin(int n, int a, vector<int> x){
		int y = 1;
		for (int i = x.size()-1; i >= 0; i--) {
			int z = y;
			y = (y*y)%n;
			if (y == 1 && z != 1 && z != (n-1)) {
				return false;
			}
			if (x[i] == 1) {
				y = (y*a)%n;
			}
		}
		if (y != 1) {
			return false;
		} else {
			return true;
		}
	}

	// performs a Miller Rabin test and stores the trace in millerRabinTrace
	bool millerRabin_verbose(int n, int a, vector<int> x){
		int y = 1;
		for (int i = x.size()-1; i >= 0; i--) {
			millerRabinTrace.push_back(x[i]);
			int z = y;
			millerRabinTrace.push_back(z);
			y = (y*y)%n;
			millerRabinTrace.push_back(y);
			millerRabinTrace.push_back(y);
			if (y == 1 && z != 1 && z != (n-1)) {
				return false;
			}
			if (x[i] == 1) {
				y = (y*a)%n;
			}
			millerRabinTrace[millerRabinTrace.size()-1] = y;
		}
		if (y != 1) {
			return false;
		} else {
			return true;
		}
	}

	void clearTrace(){
		millerRabinTrace.clear();
	}

	vector<int> getAndClearTrace(){
		vector<int> answer = millerRabinTrace;
		millerRabinTrace.clear();
		return answer;
	}

	int findPrime(){
		int answer = 0;
		bool foundPrime = false;
		int n;

		// looks for a prime number
		while (foundPrime == false) {
			n = this->getn(); // gets n
			bool nIsPossiblyPrime = true;

			// performs 20 Miller-Rabin tests
			for (int i = 1; i <=20; i++) {
				int a = this->geta(n); // gets a
				vector<int> x;
				x = this->getx(n); // gets x
	
				if (this->millerRabin(n, a, x) == false) {
					nIsPossiblyPrime = false;
					break;
				} else {
					nIsPossiblyPrime = true;
				}
			}

			if (nIsPossiblyPrime == true) {
				foundPrime = true;
			}
		}

		answer = n;
		return answer;
	}

	int gcd(int x, int y){
		int answer;
		int tempX = x;
		int tempY = y;
		x = min(tempX, tempY);
		y = max(tempX, tempY);

		int b = y % x;
		int m = (y - b) / x;

		if (b == 0) {
			answer = x;
		} else {
			answer = this->gcd(b, x);
		}

		return answer;
	}

	// performs Euclidian Algorithm
	int gcd_verbose(int x, int y){
		int answer;
		int tempX = x;
		int tempY = y;
		x = min(tempX, tempY);
		y = max(tempX, tempY);

		int b = y % x;
		int m = (y - b) / x;

		printf("%d %d %d %d \n", y, m, x, b);

		if (b == 0) {
			answer = x;
		} else {
			answer = this->gcd(b, x);
		}

		return answer;
	}

	// performs Extended Euclidian Algorithm
	int modularInverse(int e, int phi_n){
		int answer;
		if (this->gcd(e, phi_n) != 1) {
			return 0;
		}
		
		vector<int> q;
		vector<int> r;
		vector<int> s;
		vector<int> t;

		// first row
		r.push_back(phi_n);
		r.push_back(e);
		r.push_back(r[r.size()-2]%r[r.size()-1]);
		q.push_back((r[r.size()-3]-r[r.size()-1])/r[r.size()-2]);
		s.push_back(1);
		t.push_back(0);

		// second row
		r.push_back(r[r.size()-2]%r[r.size()-1]);
		q.push_back((r[r.size()-3]-r[r.size()-1])/r[r.size()-2]);
		s.push_back(0);
		t.push_back(1);

		// third row
		s.push_back(s[s.size()-2]-q[s.size()-2]*s[s.size()-1]);
		t.push_back(t[t.size()-2]-q[t.size()-2]*t[t.size()-1]);

		while (r.back() > 0) {
			r.push_back(r[r.size()-2]%r[r.size()-1]);
			q.push_back((r[r.size()-3]-r[r.size()-1])/r[r.size()-2]);
			s.push_back(s[s.size()-2]-q[s.size()-2]*s[s.size()-1]);
			t.push_back(t[t.size()-2]-q[t.size()-2]*t[t.size()-1]);
		}

		answer = this->normalize(t.back(), phi_n);

		return answer;
	}

	int modularInverse_verbose(int e, int phi_n){
		int answer;
		if (this->gcd(e, phi_n) != 1) {
			return 0;
		}
		
		vector<int> q;
		vector<int> r;
		vector<int> s;
		vector<int> t;

		// first row
		r.push_back(phi_n);
		r.push_back(e);
		r.push_back(r[r.size()-2]%r[r.size()-1]);
		q.push_back((r[r.size()-3]-r[r.size()-1])/r[r.size()-2]);
		s.push_back(1);
		t.push_back(0);

		// second row
		r.push_back(r[r.size()-2]%r[r.size()-1]);
		q.push_back((r[r.size()-3]-r[r.size()-1])/r[r.size()-2]);
		s.push_back(0);
		t.push_back(1);

		// third row
		s.push_back(s[s.size()-2]-q[s.size()-2]*s[s.size()-1]);
		t.push_back(t[t.size()-2]-q[t.size()-2]*t[t.size()-1]);

		while (r.back() > 0) {
			r.push_back(r[r.size()-2]%r[r.size()-1]);
			q.push_back((r[r.size()-3]-r[r.size()-1])/r[r.size()-2]);
			s.push_back(s[s.size()-2]-q[s.size()-2]*s[s.size()-1]);
			t.push_back(t[t.size()-2]-q[t.size()-2]*t[t.size()-1]);
		}

		for (int i = 0; i < q.size(); i++) {
			printf("%d %d %d %d %d %d\n", q[i], r[i], r[i+1], r[i+2], s[i], t[i]);
		}
		printf("    %d %d\n", s.back(), t.back());

		answer = this->normalize(t.back(), phi_n);

		return answer;
	}

	int normalize(int d, int phi_n){
		while (d <= 0) {
			d += phi_n;
		}
		while (d >= phi_n) {
			d -= phi_n;
		}
		return d;
	}

};

class CertificateGenerator{
public:
	CertificateGenerator(){}

	// performa fast exponentiation
	int fastExponentiation(int a, vector<int> x, int n){
		int y = 1;
		for (int i = 0; i < x.size(); i++) {
			y = (y*y)%n;
			if (x[i] == 1) {
				y = (a*y)%n;
			}
		}
		return y;
	}

	int fastExponentiation_verbose(int a, vector<int> x, int n){
		int y = 1;
		cout << (bitset<32>) y << endl;
		for (int i = 0; i < x.size(); i++) {
			y = (y*y)%n;
			if (x[i] == 1) {
				y = (a*y)%n;
			}
			cout << (bitset<32>) y << endl;
		}
		return y;
	}
};

class Authenticator{
private:
	int k;
public:
	Authenticator(){
		// random number seed
		srand(time(NULL));
	}

	int getu(int n){

		int answer = 0b1;

		// finds the length of u_(k-2) to u_0
		int randomBitsToGet = 0;
		n = n >> 2;
		while (n > 0) {
			n = n >> 1;
			randomBitsToGet++;
		}		

		this->k = randomBitsToGet;

		// adds the desired random bits to answer; 
		for (int i = 0; i < randomBitsToGet; i++) {
			answer = answer << 1;
			int randomInt = rand();
			randomInt = randomInt & 0b1;
			answer += randomInt;
		}

		return answer;
	}

	int getk(){
		return this->k;
	}

	// performs xor hashing algorithm
	int h(int u) {
		int answer = 0;
		for (int i = 0; i < 4; i++) {
			answer = answer ^ u;
			u = u >> 8;
		}
		answer = answer & 0x000000FF;
		return answer;
	}
};

int main(int argc, char **argv) {

	// gets n
	printf("104\n");
	RSASystem *myRSASystem = new RSASystem();
	int n = myRSASystem->getn_verbose();
	cout << "n: " << (bitset<32>) n << endl;

	// Miller-Rabin test to find non prime number
	printf("\n119\n");
	vector<int> mrTrace119;
	bool foundNonprimen = false;

	// finds a non-prime number
	while (foundNonprimen == false) {

		n = myRSASystem->getn(); // gets n

		// performs 20 Miller-Rabin tests
		for (int i = 1; i <=20; i++) {
			int a = myRSASystem->geta(n); // gets a

			vector<int> x;
			x = myRSASystem->getx(n); // gets x

			if (myRSASystem->millerRabin_verbose(n, a, x) == false) {
				vector<int> B = myRSASystem->getAndClearTrace();
				mrTrace119.insert(mrTrace119.end(), B.begin(), B.end());
				foundNonprimen = true;
				break;
			} else {
				myRSASystem->clearTrace();
			}
		}

		// prints the trace of the Miller-Rabin test
		if (foundNonprimen) {
			cout << "Found an n that is not prime." << endl;
			cout << "n: " << n << endl;
			cout << "a: " << myRSASystem->geta(n) << endl;
			printf("    x     z     y     y\n");
			for (int i = 0; i < mrTrace119.size(); i+=4) {
				printf("%5d %5d %5d %5d\n",
					mrTrace119[i],
					mrTrace119[i+1],
					mrTrace119[i+2],
					mrTrace119[i+3]);
			}
		}
		myRSASystem->clearTrace();
		mrTrace119.clear();
	}
	myRSASystem->clearTrace();

	// Miller-Rabin test to find prime number
	printf("\n123\n");
	vector<int> mrTrace123;
	bool foundPrimen = false;

	// finds a prime number
	while (foundPrimen == false) {
		n = myRSASystem->getn(); // gets n
		bool nIsPossiblyPrime = true;

		// performs 20 Miller-Rabin tests or how ever many tests until a
		// non-prime is found
		for (int i = 1; i <=20; i++) {
			int a = myRSASystem->geta(n); // gets a
			vector<int> x;
			x = myRSASystem->getx(n); // gets x

			// performs one Miller-Rabin test
			if (myRSASystem->millerRabin_verbose(n, a, x) == false) {
				myRSASystem->clearTrace();
				nIsPossiblyPrime = false;
				break;
			} else {
				nIsPossiblyPrime = true;
				mrTrace123.push_back(-1);
				mrTrace123.push_back(-1);
				mrTrace123.push_back(-1);
				mrTrace123.push_back(-1);
				vector<int> B = myRSASystem->getAndClearTrace();
				mrTrace123.insert(mrTrace123.end(), B.begin(), B.end());
			}
		}

		// checks if n is still possibly prime after 20 Miller-Rabin tests
		if (nIsPossiblyPrime == true) {
			foundPrimen = true;
		}

		// prints the trace of the 20 Miller-Rabin tests
		if (foundPrimen) {
			cout << "Found an n that is prime with high probability." << endl;
			cout << "n: " << n << endl;
			cout << "a: " << myRSASystem->geta(n) << endl;
			printf("    x     z     y     y\n");
			int testCounter = 1;
			for (int i = 0; i < mrTrace123.size(); i+=4) {
				if (mrTrace123[i] == -1) {
					cout << endl;
					cout << "Test " << testCounter << endl;
					testCounter++;
					continue;
				}
				printf("%5d %5d %5d %5d\n",
					mrTrace123[i],
					mrTrace123[i+1],
					mrTrace123[i+2],
					mrTrace123[i+3]);
			}
		}

		// clears the data structure that holds the trace
		mrTrace123.clear();
	}

	printf("\n142\n");

	// finds p and q
	int p = myRSASystem->findPrime();
	int q = myRSASystem->findPrime();
	while (p == q) {
		q = myRSASystem->findPrime();
	}

	// calculates n
	n = p * q;

	// calculates phi(n)
	int phi_n = n - 1;

	// finds e
	int e = 3;
	int gcd;
	cout << endl;
	gcd = myRSASystem->gcd_verbose(e, phi_n);
	while (gcd != 1) {
		e++;
		cout << endl;
		gcd = myRSASystem->gcd_verbose(e, phi_n);

		// repicks random primes if e < phi_n is not found
		if (e == phi_n) {
			p = myRSASystem->findPrime();
			q = myRSASystem->findPrime();
			n = p * q;
			phi_n = n - 1;
			cout << endl;
			gcd = myRSASystem->gcd_verbose(e, phi_n);
		}
	}

	// finds d
	int d = myRSASystem->modularInverse_verbose(e, phi_n);

	printf("\n156\n");
	cout << "p: " << p << endl;
	cout << "q: " << q << endl;
	cout << "n: " << n << endl;
	cout << "n: " << (bitset<32>) n << endl;
	cout << "e: " << e << endl;
	cout << "d: " << d << endl;


	printf("\n185\n");
	// creates r
	vector<uint8_t> r;

	uint8_t blank = ' '; r.push_back(blank);
	uint8_t blankA = 'A'; r.push_back(blankA);
	uint8_t blankL = 'l'; r.push_back(blankL);
	uint8_t blankI = 'i'; r.push_back(blankI);
	uint8_t blankC = 'c'; r.push_back(blankC);
	uint8_t blankE = 'e'; r.push_back(blankE);
	uint8_t n1 = n >> 24; r.push_back(n1);
	uint8_t n2 = n >> 16; r.push_back(n2);
	uint8_t n3 = n >> 8; r.push_back(n3);
	uint8_t n4 = n; r.push_back(n4);
	uint8_t e1 = e >> 24; r.push_back(e1);
	uint8_t e2 = e >> 16; r.push_back(e2);
	uint8_t e3 = e >> 8; r.push_back(e3);
	uint8_t e4 = e; r.push_back(e4);

	cout << "r: ";
	cout << (bitset<8>) blank;
	cout << (bitset<8>) blankA;
	cout << (bitset<8>) blankL;
	cout << (bitset<8>) blankI;
	cout << (bitset<8>) blankC;
	cout << (bitset<8>) blankE;
	cout << (bitset<8>) n1;
	cout << (bitset<8>) n2;
	cout << (bitset<8>) n3;
	cout << (bitset<8>) n4;
	cout << (bitset<8>) e1;
	cout << (bitset<8>) e2;
	cout << (bitset<8>) e3;
	cout << (bitset<8>) e4;
	cout << endl;

	// performs h(r)
	uint8_t h = 0;
	for (int i = 0; i < r.size(); i++) {
		h = h ^ r[i];
	}

	cout << "h(r): " << (bitset<8>) h << endl;

	CertificateGenerator *myCertificateGenerator = new CertificateGenerator();

	// vectorfies d
	vector<int> dVector;
	dVector.push_back(((d & 0x80000000) >> 31) & 0b1);
	for (int i = 30; i >= 0; i--) {
		int temp = (d & 0x7FFFFFFF) >> i;
		dVector.push_back(temp);
	}

	// calculates s
	int s = myCertificateGenerator->fastExponentiation((int) h, dVector, n);

	cout << "s: " << (bitset<32>) s << endl;
	cout << "h(r): " << (int) h << endl;
	cout << "s: " << s << endl;

	printf("\n206\n");
	Authenticator *myAuthenticator = new Authenticator();
	int u = myAuthenticator->getu(n); // gets u
	int k = myAuthenticator->getk(); // gets k
	cout << "k: " << k << endl;
	cout << "u: " << u << endl;

	printf("\n208\n");
	cout << "u: " << (bitset<32>) u << endl;

	printf("\n215\n");
	int h_u = myAuthenticator->h(u); // gets h(u)
	int v = myCertificateGenerator->fastExponentiation(h_u, dVector, n); // gets v

	// vectorfies encryption key e
	vector<int> eVector;
	eVector.push_back(((e & 0x80000000) >> 31) & 0b1);
	for (int i = 30; i >= 0; i--) {
		int temp = (e & 0x7FFFFFFF) >> i;
		eVector.push_back(temp);
	}

	// performs the encryption via fast exponentiation
	int E_RSA = myCertificateGenerator->fastExponentiation(h_u, eVector, n);

	cout << "u: " << (bitset<32>) u << endl;
	cout << "h(u): " << (bitset<8>) h_u << endl;
	cout << "v: " << (bitset<32>) v << endl;
	cout << "E_RSA(e, v): " << (bitset<32>) E_RSA << endl;

	printf("\n219\n");
	myCertificateGenerator->fastExponentiation_verbose(h_u, eVector, n);

	return 0;
}
