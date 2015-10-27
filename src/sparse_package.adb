with Ada.Text_IO, Ada.Numerics.Generic_Real_Arrays; 

package body Sparse_Package is
   package Real_Arrays is new Ada.Numerics.Generic_Real_Arrays (Real);
   

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
   function Triplet_To_Matrix (I      : in Int_Array;
			       J      : in Int_Array;
			       X      : in Real_Array;
			       N_Row  : in Pos := 0;
			       N_Col  : in Pos := 0;
			       Format : in Matrix_Format := CSC) 
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
      for K in 1 .. Int (Item'Length) loop
   	 Vector (K) := Item (K + Offset);
      end loop;
      return Vector;
   end Vectorize;
   
   function Vectorize (Item : in Int_Array) return Int_Vector is
      Vector : Int_Vector;
      Offset : constant Int := Item'First - 1;
   begin
      Vector.Set_Length (Item'Length);
      for K in 1 .. Int (Item'Length) loop
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
   function Eye (N : in Nat) return Matrix is separate;
   function Zero_Vector (N : in Nat) return Matrix is separate;
   function Dot_Product (Left_I, Right_J : in Int_Array;
			 Left_X, Right_Y : in Real_Array) return Real is separate;
   function Dot_Product_RV (X, Y : in Real_Vector) return Real is separate;

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
   function Add_RV_RV (Left, Right : in Real_Vector) return Real_Vector is separate;
   function Minus_RV_RV (Left, Right : in Real_Vector) return Real_Vector is separate;
   function Permute_By_Col (Mat : in Matrix;
			    P   : in Int_Array) return Matrix is separate;
   function Permute (Mat : in Matrix;
		     P   : in Int_Array;
		     By  : in Permute_By_Type := Column) return Matrix is separate;
   
   function Norm2_RV (X : in Real_Vector) return Real is separate;
   function Norm_RV (X : in Real_Vector) return Real is separate;
   
   function BiCGSTAB (A   : in     Matrix;
		      B   : in     Real_Vector;
		      X0  : in     Real_Vector;
		      Err :    out Real;
		      Tol : in     Real	    := 1.0e-10) 
		     return Real_Vector is separate;

   function Number_Of_Elements (X : in Matrix) return Int is (Int (X.X.Length));
   function Length (X : in Real_Vector) return Int is (Int (X.Length));
   function To_Sparse (Mat : in Matrix) return Sparse_Ptr is separate;
   
   function LU_Decomposition (Mat : in Matrix;
			      Tol : in Real   := 1.0e-12) return LU_Type is
      Sparse : Sparse_Ptr := Mat.To_Sparse;
      LU     : LU_Type;
   begin
      LU.Symbolic := CS_Sqr (Prob => Sparse);
      pragma Assert (LU.Symbolic /= null,
		     "ERROR: Problem with LU.Symbolic from CS_Sqr");
      LU.Numeric  := CS_LU (Sparse, LU.Symbolic, Tol);
      pragma Assert (LU.Numeric /= null,
		     "ERROR: Problem with LU.Numeric from CS_LU");
      LU.NCol     := Sparse.N;
      Sparse      := Free (Sparse);
      return LU;
   end LU_Decomposition;
   
   
   
   
   function Solve (LU : in LU_Type;
		   B  : in Real_Array) return Real_Array is
      X : Real_Ptrs.Pointer := Solve (LU, B);
      Y : Real_Array (B'Range) with Convention => C, Address => X.all'Address;
   begin
      return Y;
   end Solve;
   
   function Solve (LU : in LU_Type;
		   B  : in Real_Vector) return Real_Array is 
      (Solve (LU, To_Array (B)));
      
   function Solve (LU : in LU_Type;
		   B  : in Real_Vector) return Real_Vector is
      (Vectorize (Solve (LU, To_Array (B))));
      
   function Solve (LU : in LU_Type;
		   B  : in Real_Array) return Real_Ptrs.Pointer is
      use C;
      Err : C.int;
      X : Real_Ptrs.Pointer;
   begin
      X := Solve_CS (LU.NCol, LU.Symbolic, LU.Numeric, B, Err);
      pragma Assert (Err /= 1, "ERROR from Solve_CS");
      return X;
   end Solve;
      
   function Is_Valid (P	: in Real_Ptrs.Pointer;
		      N	: in Pos) return Boolean is
      X : Real_Array (1 .. N) with Convention => C, Address => P.all'Address;
   begin
      return (for all Y of X => Y'Valid);
   end Is_Valid;
   
   function Is_Valid (Mat : in Matrix) return Boolean is
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
   
   function N_Col (LU : in LU_Type) return Pos is (LU.NCol);
   
   function Max (X : in Int_Vector) return Int is
      Result : Int;
   begin
      case X.Length is
	 when 0 => return 0;
	 when 1 => return X (1);
	 when others =>
	    Result := X (1);
	    for I in 2 .. Int (X.Length) loop
	       if X (I) > Result then Result := X (I); end if;
	    end loop;
	    return Result;
      end case;
   end Max;
   
   function Max (X : in Real_Vector) return Real is
      Result : Real;
   begin
      case X.Length is
	 when 0 => return 0.0;
	 when 1 => return X (1);
	 when others =>
	    Result := X (1);
	    for I in 2 .. Int (X.Length) loop
	       if X (I) > Result then Result := X (I); end if;
	    end loop;
	    return Result;
      end case;
   end Max;
   
   function Triplet_To_Matrix (I      : in Int_Vector;
			       J      : in Int_Vector;
			       X      : in Real_Vector;
			       N_Row  : in Pos		 := 0;
			       N_Col  : in Pos		 := 0;
			       Format : in Matrix_Format := CSC) return Matrix is
      Result : Matrix;
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
				return Matrix is
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
   
   
   function Basis_Vector (I, N : in Int) return Real_Vector is
      use Ada.Containers;
      Result : Real_Vector;
   begin
      Result.Set_Length (Count_Type (N));
      Result (I) := 1.0;
      return Result;
   end Basis_Vector;
   
   procedure Print (V : in Real_Vector) is
      use Real_IO, Ada.Text_IO;
   begin
      for X of V loop
	 Put (X); New_Line;
      end loop;
   end Print;
   
   
   procedure Set_Length (V : in out Real_Vector;
			 N : in     Int) is
      use Ada.Containers;
   begin
      V.Set_Length (Count_Type (N));
   end Set_Length;
   
   
   function Abs_Max (Item : in Real_Vector) return Real is
      Result : Real := 0.0;
   begin
      for X of Item loop
	 if abs (X) > Result then Result := X; end if;
      end loop;
      return Result;
   end Abs_Max;
begin
   null;
end Sparse_Package;
