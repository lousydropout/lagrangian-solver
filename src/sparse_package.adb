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
   function N_Row (Mat : in Matrix) return Pos is separate;
   function N_Col (Mat : in Matrix) return Pos is separate;
   function Max (Item : in Int_Array) return Int is
      Result : Int := Item (Item'First);
   begin
      for N of Item loop
	 Result := Int'Max (Result, N);
      end loop;
      return Result;
   end Max;

   function Max (Item : in Real_Array) return Real is
      Result : Real := Item (Item'First);
   begin
      for N of Item loop
	 Result := Real'Max (Result, N);
      end loop;
      return Result;
   end Max;
   
   function Abs_Max (Item : in Int_Array) return Int is
      Result : Int := 0;
   begin
      for N of Item loop
	 Result := Int'Max (Result, abs (N));
      end loop;
      return Result;
   end Abs_Max;
   
   
   
   
   function Abs_Max (Item : in Real_Array) return Real is
      Result : Real := 0.0;
   begin
      for N of Item loop
	 Result := Real'Max (Result, abs (N));
      end loop;
      return Result;
   end Abs_Max;
   
   
   
   
   
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
   function Cumulative_Sum (Item : in Int_Array) return Int_Array is
      Result : Int_Array (Item'Range);
      Tmp    : Int := 1;
   begin
      for I in Item'Range loop
	 Result (I) := Tmp;
	 Tmp := Tmp + Item (I);
      end loop;
      return Result;
   end Cumulative_Sum;
   procedure Remove_Duplicates (Mat : in out Matrix) is
      N, Iter : Pos := 0;
      J : Int_Array  (1 .. Nat (Mat.P.Length)) := (others => 0);
      I : Int_Array  (1 .. Nat (Mat.I.Length));
      X : Real_Array (1 .. Nat (Mat.X.Length));
   begin
      for K in 1 .. Nat (Mat.P.Length) - 1 loop
	 Iter := 0;
	 for L in Mat.P (K) .. Mat.P (K + 1) - 1 loop
	    if Iter /= Mat.I (L) then
	       N     := N + 1;
	       Iter  := Mat.I (L);
	       I (N) := Iter; 
	       J (K) := J (K) + 1;
	       X (N) := Mat.X (L);
	    else
	       X (N) := X (N) + Mat.X (L);
	    end if;
	 end loop;
      end loop;
      J := Cumulative_Sum (J);
      Mat.I := Vectorize (I (1 .. N)); 
      Mat.X := Vectorize (X (1 .. N)); 
      Mat.P := Vectorize (J);
   end Remove_Duplicates;
   procedure Compress (Mat : in out Matrix) is
      X          : Real_Array (1 .. Int (Mat.X.Length)) := (others => 0.0);
      I          : Int_Array  (1 .. Int (Mat.X.Length)) := (others => 0);
      Col, Count : Int_Array  (1 .. Mat.N_Col + 1)          := (others => 0);
      Index      : Nat                                 := 1;
      P          : Int_Vector renames Mat.P;
   begin
      Mat.Format := CSC;
      
      for K of P loop
	 Count (K) := Count (K) + 1;
      end loop;
      
      Col := Cumulative_Sum (Count); Count := Col;

      for K in 1 .. Nat (Mat.X.Length) loop
	 Index       := Col (P (K));
	 Col (P (K)) := Col (P (K)) + 1;
	 I (Index)   := Mat.I (K); 
	 X (Index)   := Mat.X (K);
      end loop;
      
      Mat.I := Vectorize (I); 
      Mat.X := Vectorize (X);
      P     := Vectorize (Count); 
      
      Mat.Convert; Mat.Convert;
      Mat.Remove_Duplicates;
   end Compress;

   
   function Convert (Mat : in Matrix) return Matrix is
      Result : Matrix := Mat;
   begin
      return Result.Convert;
   end Convert;
   
   procedure Convert (Mat : in out Matrix) is
      X          : Real_Array (1 .. Nat (Mat.X.Length)) 
	:= (others => 0.0);
      I          : Int_Array (1 .. Nat (Mat.I.Length))  
	:= (others => 0);
      N          : constant Nat := Nat'Max (Mat.N_Col, Mat.N_Row);
      Row, Count : Int_Array (1 .. N + 1) := (others => 0);
      Index      : Int := 1;
      TRANSPOSE_EXCEPTION : exception;
   begin
      case Mat.Format is
	 when CSC => Mat.Format := CSR;
	 when CSR => Mat.Format := CSC;
	 when Triplet => raise TRANSPOSE_EXCEPTION;
      end case;
      
      for K of Mat.I loop
	 Count (K) := Count (K) + 1;
      end loop;
      
      Row := Cumulative_Sum (Count); Count := Row;
      for K in 1 .. Nat (Mat.P.Length) - 1 loop
	 for J in Mat.P (K) .. Mat.P (K + 1) - 1 loop
	    Index           := Row (Mat.I (J));
	    Row (Mat.I (J)) := Row (Mat.I (J)) + 1;
	    I (Index)       := K;
	    X (Index)       := Mat.X (J);
	 end loop;
      end loop;
      
      case Mat.Format is
	 when CSC => Mat.P := Vectorize (Count (1 .. Mat.N_Col + 1)); 
	 when CSR => Mat.P := Vectorize (Count (1 .. Mat.N_Row + 1)); 
	 when others => raise TRANSPOSE_EXCEPTION;
      end case;
      Mat.I := Vectorize (I); Mat.X := Vectorize (X);
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

   
   
   
   
begin
   null;
end Sparse_Package;
