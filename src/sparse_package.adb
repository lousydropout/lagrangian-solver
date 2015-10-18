with Ada.Text_IO, Ada.Numerics.Generic_Elementary_Functions, Ada.Numerics.Generic_Real_Arrays; 

package body Sparse_Package is
   --  package Real_Functions is
   --     new Ada.Numerics.Generic_Elementary_Functions (Real);
   package Real_Arrays is new Ada.Numerics.Generic_Real_Arrays (Real);
   package Int_IO is new Ada.Text_IO.Integer_IO (Integer);
   package Real_IO is new Ada.Text_IO.Float_IO (Real);
   package Matrix_Format_IO is new Ada.Text_IO.Enumeration_IO (Matrix_Format);
   

   procedure Print (Mat : in Matrix) is separate;
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Basic Getter Functions -----------------------------------
   function Norm2 (Item : in Matrix) return Real is separate;
   function N_Row (Mat : in Matrix) return Pos is separate;
   function N_Col (Mat : in Matrix) return Pos is separate;
   function Max_Int_Array (Item : in Int_Array) return Int is separate;
   function Max_Real_Array (Item : in Real_Array) return Real is separate;
   function Abs_Max_IA (Item : in Int_Array) return Int is separate;
   function Abs_Max_RA (Item : in Real_Array) return Real is separate;
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Functions for Creating Sparse Matrices -------------------
   function Triplet_To_Matrix (I      : Int_Array;
			       J      : Int_Array;
			       X      : Real_Array;
			       N_Row  : Pos := 0;
			       N_Col  : Pos := 0;
			       Format : Matrix_Format := CSC) 
			      return Matrix is separate;
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   -------- Essential Tools -----------------------------------------
   function Cumulative_Sum (Item : in Int_Array) return Int_Array is separate;
   procedure Remove_Duplicates (Mat : in out Matrix) is separate;
   procedure Compress (Mat : in out Matrix) is separate;
   procedure Convert (Mat : in out Matrix) is separate;
   function Convert (Mat : in Matrix) return Matrix is
      Result : Matrix := Mat;
   begin
      Result.Convert;
      return Result;
   end Convert;
   
   
   -- Vectorize & To_Array are needed in Triplet_To_Matrix
   function Vectorize (Item : in Real_Array) return Real_Vector is
      Vector : Real_Vector;
      Offset : constant Int := Item'First - 1;
   begin
      Vector.Set_Length (Item'Length);
      for K in 1 .. Item'Length loop
   	 Vector (K) := Item (K + Offset);
      end loop;
      return Vector;
   end Vectorize;
   
   function Vectorize (Item : in Int_Array) return Int_Vector is
      Vector : Int_Vector;
      Offset : constant Int := Item'First - 1;
   begin
      Vector.Set_Length (Item'Length);
      for K in 1 .. Item'Length loop
   	 Vector (K) := Item (K + Offset);
      end loop;
      return Vector;
   end Vectorize;
   
   function To_Array (Item : in Real_Vector) return Real_Array is
      Result : Real_Array (1 .. Nat (Item.Length));
   begin
      for K in Result'Range loop
	 Result (K) := Item (K);
      end loop;
      return Result;
   end To_Array;
   
   function To_Array (Item : in Int_Vector) return Int_Array is
      Result : Int_Array (1 .. Nat (Item.Length));
   begin
      for K in Result'Range loop
	 Result (K) := Item (K);
      end loop;
      return Result;
   end To_Array;

   
   function Vectorize (I : in Int_Array;
		       X : in Real_Array) return Matrix is
      Result   : Matrix;
      Offset_I : constant Int := I'First - 1;
      Offset_X : constant Int := X'First - 1;
   begin
      Result.Format := CSC;
      Result.N_Row := I (I'Last);
      Result.N_Col := 1;
      Result.P.Set_Length (2); 
      Result.I.Set_Length (I'Length);
      Result.X.Set_Length (X'Length);
      
      Result.P (1) := 1; 
      Result.P (2) := Nat (X'Length) + 1;
      for K in I'Range loop
	 Result.I (K) := I (K + Offset_I);
	 Result.X (K) := X (K + Offset_X);
      end loop;
      return Result;
   end Vectorize;

   
   
   
   
   
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Testing Functions -----------------------------------
   function Is_Col_Vector (A : in Matrix) return Boolean is separate;
   function Is_Square_Matrix (A : in Matrix) return Boolean is separate;
   function Has_Same_Dimensions (Left, Right : in Matrix) return Boolean is separate;   
   
   
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Matrix Operations -----------------------------------
   function Eye (N : in Pos) return Matrix is separate;
   function Zero_Vector (N : in Nat) return Matrix is separate;
   function Dot_Product (Left_I, Right_J : in Int_Array;
		      Left_X, Right_Y : in Real_Array) return Real is separate;
   procedure Transposed (Mat : in out Matrix) is separate;
   function Transpose (Mat : in Matrix) return Matrix is separate;
   function Mult (Left, Right : in Matrix) return Matrix is separate;
   function Mult_Int_Array (Left, Right : in Int_Array) return Boolean is separate;
   function Plus (Left  : in Matrix;
		  Right : in Matrix) return Matrix is separate;
   function Minus (Left  : in Matrix;
		   Right : in Matrix) return Matrix is separate;
   function Kronecker (Left, Right : in Matrix) return Matrix is separate;
   function Direct_Sum (Left, Right : in Matrix) return Matrix is separate;
   function Mult_R_RV (Left  : in Real;
		 Right : in Real_Vector) return Real_Vector is separate;
   function Mult_M_RV (Left  : in Matrix;
		       Right : in Real_Vector) return Real_Vector is separate;

begin
   null;
end Sparse_Package;
