%to find a target x (beta) value based on two x,y coordinate pairs and a target y (success probability) value 


    x1 = spm_input('slope: enter x1',2, 'r');
    y1 = spm_input('slope: enter y1',2, 'r');
    x2 = spm_input('slope: enter x2',2, 'r');
    y2 = spm_input('slope: enter y2',2, 'r');
    tprob = spm_input('target prob.:',2,'r');
    
    coord = [x1 y1 x2 y2]
    slope = (y2-y1)/(x2-x1)
    b = y1 - x1*slope;
    fprintf(['probability (y) = <target> when <x> has a value of: \n']);
    xval = (tprob - b)/slope

    