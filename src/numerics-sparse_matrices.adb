package body Numerics.Sparse_Matrices is
   
   procedure Print (Mat : in Sparse_Matrix) is separate;
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Basic Getter Functions -----------------------------------
   function Norm2 (Item : in Sparse_Matrix) return Real is separate;
   function N_Row (Mat : in Sparse_Matrix) return Pos is separate;
   function N_Col (Mat : in Sparse_Matrix) return Pos is separate;
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Functions for Creating Sparse Matrices -------------------
   function Sparse (X : in Real_Matrix) return Sparse_Matrix is
      use Ada.Containers;
      Result : Sparse_Matrix;
      N : constant Int := X'Length (1) * X'Length (2);
   begin
      Result.N_Row := X'Length (1);
      Result.N_Col := X'Length (2);
      Result.Format := Triplet;
      Result.I.Reserve_Capacity (Count_Type (N));
      Result.P.Reserve_Capacity (Count_Type (N));
      Result.X.Reserve_Capacity (Count_Type (N));
      for I in 1 .. Int (X'Length (1)) loop
	 for J in 1 .. Int (X'Length (2)) loop
	    Result.I.Append (I); 
	    Result.P.Append (J);
	    Result.X.Append (X (I, J));
	 end loop;
      end loop;
      Compress (Result);
      return Result;
   end Sparse;
      
   function Triplet_To_Matrix (I      : in Int_Array;
			       J      : in Int_Array;
			       X      : in Real_Array;
			       N_Row  : in Pos := 0;
			       N_Col  : in Pos := 0;
			       Format : in Sparse_Matrix_Format := CSC) 
			      return Sparse_Matrix is separate;
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   -------- Essential Tools -----------------------------------------
   --  function Cumulative_Sum (Item : in Int_Array) return Int_Array is separate;
   procedure Remove_Duplicates (Mat : in out Sparse_Matrix) is separate;
   procedure Compress (Mat : in out Sparse_Matrix) is separate;
   procedure Convert (Mat : in out Sparse_Matrix) is separate;
   function Convert (Mat : in Sparse_Matrix) return Sparse_Matrix is
      Result : Sparse_Matrix := Mat;
   begin
      Result.Convert;
      return Result;
   end Convert;
   
   
   
   function Vectorize (I : in Int_Array;
		       X : in Real_Array) return Sparse_Matrix is
      Result   : Sparse_Matrix;
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
   function Is_Col_Vector (A : in Sparse_Matrix) return Boolean is separate;
   function Is_Square_Matrix (A : in Sparse_Matrix) return Boolean is separate;
   function Has_Same_Dimensions (Left, Right : in Sparse_Matrix) return Boolean is separate;   
   
   
   
   
   ------------------------------------------------------------------
   ------------------------------------------------------------------
   ------- Matrix Operations -----------------------------------
   function Eye (N : in Nat) return Sparse_Matrix is separate;
   function Zero_Vector (N : in Nat) return Sparse_Matrix is separate;

   procedure Transposed (Mat : in out Sparse_Matrix) is separate;
   function Transpose (Mat : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Mult (Left, Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Plus (Left  : in Sparse_Matrix;
		  Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Minus (Left  : in Sparse_Matrix;
		   Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Kronecker (Left, Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Direct_Sum (Left, Right : in Sparse_Matrix) return Sparse_Matrix is separate;
   function Mult_M_RV (Left  : in Sparse_Matrix;
		       Right : in Real_Vector) return Real_Vector is separate;
   function Permute_By_Col (Mat : in Sparse_Matrix;
			    P   : in Int_Array) return Sparse_Matrix is separate;
   function Permute (Mat : in Sparse_Matrix;
		     P   : in Int_Array;
		     By  : in Permute_By_Type := Column) return Sparse_Matrix is separate;
   
   
   function BiCGSTAB (A   : in     Sparse_Matrix;
		      B   : in     Real_Vector;
		      X0  : in     Real_Vector;
		      Err :    out Real;
		      Tol : in     Real	    := 1.0e-10) 
		     return Real_Vector is separate;

   function Number_Of_Elements (X : in Sparse_Matrix) return Int is (Int (X.X.Length));
   
   function Is_Valid (Mat : in Sparse_Matrix) return Boolean is
      use IV_Package, RV_Package;
   begin
      if Mat.I = IV_Package.Empty_Vector 
	or else Mat.P = IV_Package.Empty_Vector
	or else Mat.X = RV_Package.Empty_Vector 
	or else Mat.N_Row = 0 or else Mat.N_Col = 0 
      then
	 return False;
      end if;
      return True;
   end Is_Valid;
   
   function Triplet_To_Matrix (I      : in Int_Vector;
			       J      : in Int_Vector;
			       X      : in Real_Vector;
			       N_Row  : in Pos		 := 0;
			       N_Col  : in Pos		 := 0;
			       Format : in Sparse_Matrix_Format := CSC) 
			      return Sparse_Matrix is
      Result : Sparse_Matrix;
   begin
      Result.N_Row  := (if N_Row = 0 then Max (I) else N_Row);
      Result.N_Col  := (if N_Col = 0 then Max (J) else N_Col);
      
      Result.Format := Triplet;
      Result.I := I; Result.P := J; Result.X := X;
      case Format is
	 when CSC     => 
	    Result.Compress;
	 when CSR     => 
	    Result.Compress; 
	    Result.Convert;
	 when Triplet => 
	    null;
      end case;
      return Result;
   end Triplet_To_Matrix;
   
   function Read_Sparse_Triplet (File_Name : in String;
				 Offset	   : in Int    := 0)
				return Sparse_Matrix is
      use Ada.Text_IO, Ada.Containers, Real_IO, Int_IO;
      N_Lines : Int := 0;
      I_Vec : Int_Vector;
      J_Vec : Int_Vector;
      X_Vec : Real_Vector;
      Int_Input : Int;
      Real_Input : Real;
      File : File_Type;
   begin
      Open (File => File, Mode => In_File, Name => File_Name);
      --- Count number of lines
      while not End_Of_File (File) loop 
	 Skip_Line (File); N_Lines := N_Lines + 1;
      end loop;
      Reset (File); -- Jump back to beginning of input file
      
      -- Set lengths of vectors
      I_Vec.Set_Length (Count_Type (N_Lines));
      J_Vec.Set_Length (Count_Type (N_Lines));
      X_Vec.Set_Length (Count_Type (N_Lines));
      
      for K in 1 .. N_Lines loop
	 Get (File, Int_Input); I_Vec (K)  := Int_Input + 1 - Offset;
	 Get (File, Int_Input); J_Vec (K)  := Int_Input + 1 - Offset;
	 Get (File, Real_Input); X_Vec (K) := Real_Input;
      end loop;
      Close (File);
      
      return Triplet_To_Matrix (I_Vec, J_Vec, X_Vec);
   end Read_Sparse_Triplet;
   
   
   procedure Cumulative_Sum (Item : in out Int_Array) is separate;
   
end Numerics.Sparse_Matrices;
