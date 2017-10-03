with Numerics, Numerics.Dense_Matrices;
use  Numerics, Numerics.Dense_Matrices;
generic
   Num : in Nat;
package Dense_AD is
   

   type Evaluation_Level is (Value, Gradient, Hessian);
   
   subtype Vector is Real_Vector (1 .. Num);
   subtype Matrix is Real_Matrix (1 .. Num, 1 .. Num);
   type AD_Type is private;
   type AD_Vector is array (1 .. Num) of AD_Type;
   
   procedure Set_Evaluation_Level (Value : in Evaluation_Level)
     with Inline => True;
   function Get_Evaluation_Level return Evaluation_Level with Inline => True;
   
   function Var (X  : in Real;
		 I  : in Nat;
		 Dx : in Real := 1.0) return AD_Type with Pre => I <= Num;
   
   function Const (X : in Real) return AD_Type;
   function Var (X : in Vector) return AD_Vector
     with Pre => X'Length = Num;
   
   function Val (X : in AD_Type) return Real;
   function Grad (X : in AD_Type) return Vector;
   function Hessian (X : in AD_Type) return Matrix;
   
   function "+" (X, Y : in AD_Type) return AD_Type;
   function "-" (X, Y : in AD_Type) return AD_Type;
   function "*" (X, Y : in AD_Type) return AD_Type;
   function "/" (X, Y : in AD_Type) return AD_Type;
   function "**" (X : in AD_Type; K : in Integer) return AD_Type;
   function "**" (X : in AD_Type; K : in Real) return AD_Type;
   function Sqrt (X : in AD_Type) return AD_Type is (X ** 0.5);
   function Square (X : in AD_Type) return AD_Type is (X ** 2);
   function Cube (X : in AD_Type) return AD_Type is (X ** 3);
   
   function "+" (X : in AD_Type) return AD_Type is (X) with Inline => True;
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
    
   function Sin (X : in AD_Type) return AD_Type;
   function Cos (X : in AD_Type) return AD_Type;
   function Tan (X : in AD_Type) return AD_Type;
   function Sec (X : in AD_Type) return AD_Type is (1.0 / Sin (X));
   function Csc (X : in AD_Type) return AD_Type is (1.0 / Cos (X));
   function Cot (X : in AD_Type) return AD_Type is (1.0 / Tan (X));
   
   function Exp  (X : in AD_Type) return AD_Type;
   function Log  (X : in AD_Type) return AD_Type; 
   function Sinh (X : in AD_Type) return AD_Type;
   function Cosh (X : in AD_Type) return AD_Type;
   function Tanh (X : in AD_Type) return AD_Type;
   function Sech (X : in AD_Type) return AD_Type is (1.0 / Sinh (X));
   function Csch (X : in AD_Type) return AD_Type is (1.0 / Cosh (X));
   function Coth (X : in AD_Type) return AD_Type is (1.0 / Tanh (X));
   
   function Sign (X : in AD_Type) return Real;
   
   ------------- procedures ----------------------
   procedure Print (X : in AD_Type);
   
private
   
   Level : Evaluation_Level := Hessian;
   
   type AD_Type is
      record
	 Val : Real;
	 Grad : Vector;
	 Hessian : Matrix;
      end record;
   
   G0 : constant Vector := (others => 0.0);
   H0 : constant Matrix := (others => (others => 0.0));
end Dense_AD;
