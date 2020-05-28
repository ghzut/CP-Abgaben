class Complex
{
private:
   double m_real;
   double m_imag;

public:
   // c'tors
   Complex ();
   Complex (double, double);

   // conversion constructor
   Complex (double);

   // setter / getter methods
   void SetReal (double real) { m_real = real; }
   double GetReal () { return m_real; }
   void SetImag (double imag) { m_imag = imag; }
   double GetImag () { return m_imag; }

   // unary arithmetic operators
   friend Complex operator- (const Complex&);

   // binary arithmetic operators
   friend Complex operator+ (const Complex&, const Complex&);
   friend Complex operator- (const Complex&, const Complex&);
   friend Complex operator* (const Complex&, const Complex&);
   friend Complex operator/ (const Complex&, const Complex&);

   // comparison operators
  friend bool operator== (const Complex&, const Complex&);
  friend bool operator!= (const Complex&, const Complex&);

  // complex conjunction operator
  friend Complex operator! (const Complex&);

  // input / output operators
  friend ostream& operator<< (ostream&, const Complex&);
  friend istream& operator>> (istream&, Complex&);
};
