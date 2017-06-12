with Numerics, Numerics.Sparse_Matrices;
use  Numerics, Numerics.Sparse_Matrices;

package Forward_AD.AD2D is
   type AD2D is tagged private;
   type AD2D_Vector is array (Nat range <>) of AD2D;
   
   function "+" (A, B : in AD2D) return AD2D;
   function "-" (A, B : in AD2D) return AD2D;
   function "*" (A : in AD2D;
		 B : in AD2D) return AD_Type;
   
   function "-" (A : in AD2D) return AD2D;
   function Square (A : in AD2D) return AD_Type is (A * A);
   
   
   function To_AD2D_Vector (Pos_Vector : in Pos2D_Vector) return AD2D_Vector
     with Pre => Pos_Vector'Length mod 2 = 0;
   
private
   
   type AD2D is tagged
      record
	 X : AD_Type;
	 Y : AD_Type;
      end record;
   
end Forward_AD.AD2D;
