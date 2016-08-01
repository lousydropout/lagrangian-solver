with Ada.Text_IO;
package body Forward_AD is
   
   procedure Print (X : in AD_Type) is
      use Ada.Text_IO;
   begin
      Put_Line ("   value: " & Real'Image (X.Val));
      Print (X.Grad);
      New_Line;
   end Print;
      
   
   function Var (X    : in Real;
		 I, N : in Nat;
		 Dx   : in Real	:= 1.0) return AD_Type is
      Result : AD_Type;
   begin
      Result.N   := N;
      Result.Val := X;
      Set_Length (Result.Grad, N);
      Set (Result.Grad, I, Dx);
      return Result;
   end Var;
   
   function Zero (N : in Nat) return AD_Type is
      Result : AD_Type;
   begin
      Result.N := N;
      Result.Val := 0.0;
      Set_Length (Result.Grad, N);
      return Result;
   end Zero;
   
   
   function Var (X : in Real_Array) return AD_Vector is
      Result : AD_Vector (1 .. X'Length);
   begin
      for I in X'Range loop
	 Result (I) := Var (X (I), I - X'First + 1, X'Length);
      end loop;
      return Result;
   end Var;
   
   function Val (X : in AD_Type) return Real is (X.Val);
   function Grad (X : in AD_Type) return Sparse_Vector is (X.Grad);
   function Grad (X : in AD_Type) return Real_Array is (To_Array (X.Grad));
   function Length (X : in AD_Type) return Pos is (X.N);
   function "+" (X, Y : in AD_Type) return AD_Type is
   begin
      return (N    => X.N, 
	      Val  => X.Val  + Y.Val, 
	      Grad => X.Grad + Y.Grad);
   end "+";
   function "-" (X, Y : in AD_Type) return AD_Type is
   begin
      return (N    => X.N, 
	      Val  => X.Val  - Y.Val, 
	      Grad => X.Grad - Y.Grad);
   end "-";
   function "*" (X, Y : in AD_Type) return AD_Type is
   begin
      return (N    => X.N,
	      Val  => X.Val * Y.Val,
	      Grad => X.Val * Y.Grad + Y.Val * X.Grad);
   end "*";
   
   function "*" (X : in Real; Y : in AD_Vector) return AD_Vector is
      Z : AD_Vector := Y;
   begin
      for W of Z loop
	 W := X * W;
      end loop;
      return Z;
   end "*";
   
   function "*" (X : in Real_Matrix; Y : in AD_Vector) return AD_Vector is
      Z : AD_Vector := 0.0 * Y;
   begin
      for K in Y'Range loop
	 for I in X'Range (1) loop
	    Z (I) := Z (I) + X (I, K) * Y (K);
	 end loop;
      end loop;
      return Z;
   end "*";
   
   
   function "/" (X, Y : in AD_Type) return AD_Type is
      Z : constant Real := 1.0 / Y.Val;
   begin
      return (N    => X.N,
	      Val  => X.Val * Z,
	      Grad => Z * X.Grad - (Z ** 2 * X.Val) * Y.Grad);
   end "/";
   
   function "**" (X : in AD_Type; N : Pos) return AD_Type is
   begin
      return (N    => X.N,
   	      Val  => X.Val ** Natural (N),
   	      Grad => (Real (N) * X.Val ** Natural (N - 1)) * X.Grad);
   end "**";
   function Sin (X : in AD_Type) return AD_Type is
      use Real_Functions;
   begin
      return (N    => X.N, 
	      Val  => Sin (X.Val), 
	      Grad => Cos (X.Val) * X.Grad);
   end Sin;
   function Cos (X : in AD_Type) return AD_Type is
      use Real_Functions;
   begin
      return (N    =>  X.N, 
	      Val  =>  Cos (X.Val), 
	      Grad => -Sin (X.Val) * X.Grad);
   end Cos;
   function Tan (X : in AD_Type) return AD_Type is
      use Real_Functions;
      Y     : constant Real := Cos (X.Val);
      T_Val : Real;
   begin
      pragma Assert (Y /= 0.0);
      T_Val := Sin (X.Val) / Y;
      return (N    => X.N,
	      Val  => T_Val,
	      Grad => (1.0 - T_Val ** 2) * X.Grad);
   end Tan;
	      
      
   function Exp (X : in AD_Type) return AD_Type is
      use Real_Functions;
      Y : constant Real := Exp (X.Val);
   begin
      return (N    => X.N,
	      Val  => Y,
	      Grad => Y * X.Grad);
   end Exp;
   
   
   function Log (X : in AD_Type) return AD_Type is
      use Real_Functions;
   begin
      pragma Assert (X.Val > 0.0);
      return (N    => X.N,
	      Val  => Log (X.Val),
	      Grad => (1.0 / X.Val) * X.Grad);
   end Log;
   
   
   
   -----------------------------------------------------------
   function "+" (X : in AD_Type) return AD_Type is
   begin
      return X;
   end "+";
   function "-" (X : in AD_Type) return AD_Type is
   begin
      return (N    =>  X.N, 
	      Val  => -X.Val,
	      Grad => -X.Grad);
   end "-";
   
   
   
   
   function "*" (X : in AD_Type; Y : in Real) return AD_Type is
   begin
      return (N    => X.N,
	      Val  => Y * X.Val,
	      Grad => Y * X.Grad);
   end "*";
   function "*" (Y : in Real; X : in AD_Type) return AD_Type is
   begin
      return (N    => X.N,
	      Val  => Y * X.Val,
	      Grad => Y * X.Grad);
   end "*";
   function "/" (X : in AD_Type; Y : in Real) return AD_Type is
   begin
      return (1.0 / Y) * X;
   end "/";

   
end Forward_AD;
