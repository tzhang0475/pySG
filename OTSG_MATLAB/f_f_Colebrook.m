function f = f_f_Colebrook(Re)
ff=@(f) (1/sqrt(f)+2.0*log10(2.51/(Re*sqrt(f))));
if Re <= 2300
    f = 64 / Re;
else if Re <= 1.0e4
        f = fzero(ff,0.036);
    else if Re <= 1.0e5
            f = fzero(ff,0.024);
        else if Re <= 1.0e6
                f = fzero(ff,0.015);
            else if Re <= 1.0e7
                    f = fzero(ff,0.009);
                else
                    f = fzero(ff,0.005);
                end
            end
        end
    end
end