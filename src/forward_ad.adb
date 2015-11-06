package body Forward_AD is
   
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
   
   function Val (X : in AD_Type) return Real is (X.Val);
   function Grad (X : in AD_Type) return Sparse_Vector is (X.Grad);
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
   function Exp (X : in AD_Type) return AD_Type is
      use Real_Functions;
      Y : constant Real := Exp (X.Val);
   begin
      return (N    => X.N,
	      Val  => Y,
	      Grad => Y * X.Grad);
   end Exp;
end Forward_AD;
