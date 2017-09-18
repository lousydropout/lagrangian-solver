with Numerics, Numerics.Dense_Matrices;
use  Numerics, Numerics.Dense_Matrices;
generic
   N : in Nat;
package Dense_AD is
   type Evaluation_Level is (Value, Gradient, Hessian);
   Level : Evaluation_Level := Hessian;
   
   subtype Vector is Real_Vector (1 .. N);
   subtype Matrix is Real_Matrix (1 .. N, 1 .. N);
   type AD_Type is private;
   type AD_Vector is array (1 .. N) of AD_Type;
   
   function Var (X  : in Real;
		 I  : in Nat;
		 Dx : in Real := 1.0) return AD_Type with Pre => I <= N;
   
   function Const (X : in Real) return AD_Type;
   function Var (X : in Vector) return AD_Vector
     with Pre => X'Length = N;
   
   function Val (X : in AD_Type) return Real;
   function Grad (X : in AD_Type) return Vector;
   function Hessian (X : in AD_Type) return Matrix;
   
   function "+" (X, Y : in AD_Type) return AD_Type;
   function "-" (X, Y : in AD_Type) return AD_Type;
   function "*" (X, Y : in AD_Type) return AD_Type;
   function "/" (X, Y : in AD_Type) return AD_Type;
   function "**" (X : in AD_Type; K : in Integer) return AD_Type;
   
   function Sin (X : in AD_Type) return AD_Type;
   function Cos (X : in AD_Type) return AD_Type;
   function Tan (X : in AD_Type) return AD_Type;
   function Exp (X : in AD_Type) return AD_Type;
   function Log (X : in AD_Type) return AD_Type; 
   
   function "+" (X : in AD_Type) return AD_Type is (X);
   function "-" (X : in AD_Type) return AD_Type;
   
   function "+" (X : in Real; Y : in AD_Type) return AD_Type;
   function "+" (X : in AD_Type; Y : in Real) return AD_Type is (Y + X);
   
   function "-" (X : in Real; Y : in AD_Type) return AD_Type;
   function "-" (X : in AD_Type; Y : in Real) return AD_Type;
   
   function "*" (X : in Real; Y : in AD_Type) return AD_Type;
   function "*" (X : in AD_Type; Y : in Real) return AD_Type is (Y * X);
   
   function "/" (X : in Real; Y : in AD_Type) return AD_Type is (X * (Y ** (-1)));
   function "/" (X : in AD_Type; Y : in Real) return AD_Type is ((1.0 / Y) * X)
     with Pre => Y /= 0.0;
		 
   ------------- procedures ----------------------
   procedure Print (X : in AD_Type);
   
private
   
   type AD_Type is
      record
	 Val : Real;
	 Grad : Vector;
	 Hessian : Matrix;
      end record;
   
   G0 : constant Vector := (others => 0.0);
   H0 : constant Matrix := (others => (others => 0.0));
end Dense_AD;
