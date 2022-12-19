classdef Noise
    %% line1
    %  line2
    %  
    %  Created 28-Sep-2022 19:37:14 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 9.13.0.2049777 (R2022b) for MACI64.  Copyright 2022 John J. Lee.

    % N.B.:  Julius O. Smith, Spectral Autio Signal Processing, W3K Publishing.
    %
    % Printed book version:
    % Smith, Julius O. 
    % Spectral Audio Signal Processing,
    % W3K Publishing, http://books.w3k.org/,
    % ISBN 978-0-9745607-3-1.
    %
    % Web version:
    % Smith, J.O. Spectral Audio Signal Processing,
    % http://ccrma.stanford.edu/~jos/sasp/, online book,
    % 2011 edition,
    % accessed <date>.
    %
    % Specific page citation example:
    % Smith, J.O. "Hamming Window", in 
    % Spectral Audio Signal Processing,
    % http://ccrma.stanford.edu/~jos/sasp/Hamming_Window.html, online book, 
    % 2011 edition,
    % accessed <date>.
    %
    % Raw HTML example:
    % Smith, J.O. 
    % Spectral Audio Signal Processing,
    % <A HREF="http://ccrma.stanford.edu/~jos/sasp/">
    % <tt>http://ccrma.stanford.edu/~jos/sasp/</tt></A>, 
    % online book, 2011 edition,
    % accessed <date>.
    %
    % BibTeX example (requires \usepackage{html} where html.sty comes from the latex2html distribution):
    % @BOOK{SASPWEB2011,
    % 	AUTHOR = "Julius O. Smith",
    % 	TITLE = "Spectral Audio Signal Processing",
    % 	PUBLISHER = "\htmladdnormallink{\texttt{http:}}{http://ccrma.stanford.edu/~jos/sasp/}\texttt{//\-ccrma.stanford.edu/\-\~{}jos/\-sasp/}",
    % 	YEAR = "accessed <date>",
    %         NOTE = "online book, 2011 edition"
    % }
    %
    % LaTeX citation example:
    % I like to cite the online book and add a footnote to the specific page, e.g.,
    % 
    % \cite{SASP}\footnote{\texttt{http://ccrma.stanford.edu/\~{}jos/sasp/Hamming\_Window.html}}
    % Or, if you want live links in the HTML version of your own online material,
    % \cite{SASP}\footnote{\htmladdnormallink{\texttt{%
    % http://ccrma.stanford.edu/\~{}jos/sasp/Hamming\_Window.html}}{%
    % http://ccrma.stanford.edu/\~{}jos/sasp/Hamming\_Window.html}}    
    
    methods
        function this = Noise(varargin)
            %% NOISE 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            ip = inputParser;
            addParameter(ip, "arg1", [], @(x) true)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
        end
    end

    methods (Static)
        function x = pink(Nx)
            %% Example: Synthesis of 1/F Noise (Pink Noise), 
            %  Pink noise or ``1/f noise'' is an interesting case because it occurs often in nature [294],7.11
            %  is often preferred by composers of computer music, and there is no exact (rational, finite-order) 
            %  filter which can produce it from white noise. This is because the ideal amplitude response of the 
            %  filter must be proportional to the irrational function $ 1/\sqrt{f}$ , where $ f$ denotes 
            %  frequency in Hz. However, it is easy enough to generate pink noise to any desired degree of 
            %  approximation, including perceptually exact.
            %  See also:  https://ccrma.stanford.edu/~jos/sasp/Example_Synthesis_1_F_Noise.html

            arguments
                Nx {mustBeScalarOrEmpty} = 2^16 % number of samples to synthesize
            end
             
            % The following Matlab/Octave code generates pretty good pink noise:
            B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
            A = [1 -2.494956002   2.017265875  -0.522189400];
            nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
            v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
            x = filter(B,A,v);    % Apply 1/F roll-off to PSD
            x = x(nT60+1:end);    % Skip transient response            
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
