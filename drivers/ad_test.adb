with Forward_AD, Ada.Text_IO, Numerics;
use  Forward_AD, Ada.Text_IO, Numerics;

procedure AD_Test is
   A, B : AD_Type;

begin
   
   A := Var (30.0, 1, 2);
   A := Exp (+A);
   --  B := A * A; -- Var (3.0, 2, 2);

   --  Put_Line ("------ A ----------------------");
   Print (A);
   --  Put_Line ("------ B ----------------------");
   --  Print (B);
   Put_Line ("------ A * B ----------------------");
   Print (A);
   null;
end AD_Test;
