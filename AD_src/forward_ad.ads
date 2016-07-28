with Numerics, Numerics.Sparse_Matrices;
use  Numerics, Numerics.Sparse_Matrices;

package Forward_AD is
   type AD_Type is tagged private;
     
   function Var (X    : in Real;
		 I, N : in Nat;
		 Dx   : in Real	:= 1.0) return AD_Type;
   
   function Val (X : in AD_Type) return Real;
   function Grad (X : in AD_Type) return Sparse_Vector;
   function Length (X : in AD_Type) return Pos;
   function "+" (X, Y : in AD_Type) return AD_Type with Pre => X.Length = Y.Length;
   function "-" (X, Y : in AD_Type) return AD_Type with Pre => X.Length = Y.Length;
   function "*" (X, Y : in AD_Type) return AD_Type with Pre => X.Length = Y.Length;
   function Sin (X : in AD_Type) return AD_Type;
   function Cos (X : in AD_Type) return AD_Type;
   function Exp (X : in AD_Type) return AD_Type;
   
private
   type AD_Type is tagged
      record
	 N    : Pos := 0;
	 Val  : Real;
	 Grad : Numerics.Sparse_Vector;
      end record;
end Forward_AD;
