function result=clarify_design_horizon()
%CLARIFY_DESIGN_HORIZON Post-run clarification, preserving the 30 s results.
root=fileparts(fileparts(mfilename('fullpath')));out=fullfile(root,'results','design_sensitivity');
T=readtable(fullfile(out,'all_plant_cases.csv'),'TextType','string');D=readtable(fullfile(out,'design_register.csv'),'TextType','string');
rows=T(T.Stable & T.HorizonLimited,:);assert(height(rows)==1,'Reassess scope if horizon-limited cases change.');r=rows(1,:);c=D(D.DesignID==r.DesignID,:);
source=jsondecode(fileread(fullfile(root,'results','research_summary.json')));p=source.p;
names={'Tg','Tt','Tp','Ksys','R'};for j=1:5,p.(names{j})=p.(names{j})*r.([names{j} '_factor']);end
assert(c.Type=="PI");A=[-1/p.Tg 0 -1/(p.R*p.Tg);1/p.Tt -1/p.Tt 0;0 p.Ksys/p.Tp -1/p.Tp];B=[1/p.Tg;0;0];E=[0;0;-p.Ksys/p.Tp];C=[0 0 1];
Ac=[A-B*c.Kp*C B*c.Ki;-C 0];Ec=[E;0];Cy=[C 0];Cu=[-c.Kp*C c.Ki];
[V,lam]=eig(Ac,'vector');xss=-Ac\(Ec*p.loadStep);res=-(Cy*V).'.*(V\xss);fss=Cy*xss;
freq=@(tau)real(exp(tau(:)*lam.')*res+fss);band=source.specifications.recoveryBandHz;
tau=(0:.001:60)';y=freq(tau);ys=lsim(ss(Ac,Ec,Cy,0),p.loadStep*ones(size(tau)),tau);err=max(abs(y-ys));assert(err<1e-10);
fine=(0:.0001:60)';yf=freq(fine);last=find(abs(y)>band,1,'last');lastf=find(abs(yf)>band,1,'last');assert(last<numel(tau) && lastf<numel(fine));
fun=@(t)abs(freq(t))-band;cross=fzero(fun,[tau(last),tau(last+1)]);crossFine=fzero(fun,[fine(lastf),fine(lastf+1)]);
tail=abs(fss)+sum(abs(res).*exp(real(lam)*60));assert(tail<band && abs(cross-crossFine)<1e-8);
result=struct('scope','Post-run numerical clarification of the single horizon-limited row; original 30s primary results unchanged.', ...
 'DesignID',char(r.DesignID),'CaseID',r.CaseID,'originalElapsedEndpoint_s',29,'originalFrequencyAtEndpoint_Hz',freq(29), ...
 'originalRecovery_s',r.Recovery_s,'extendedElapsedEndpoint_s',60,'extendedSampleInterval_s',.001, ...
 'extendedSampledRecovery_s',tau(last+1),'refinedRecovery_s',cross,'fineGridRootDifference_s',crossFine-cross, ...
 'modalLsimMaxDifference_Hz',err,'tailBoundAtExtendedEndpoint_Hz',tail,'originalAllTargetsPass',logical(r.AllTargetsPass), ...
 'allTargetsClassificationChanged',false,'eligibleMatchedRankingsChanged',false, ...
 'reason','Recovery remains far above 8 s; this unpaired PI design has no positive-gain PIDF counterpart under the fixed closure.');
writetable(table(tau,y,ys,'VariableNames',{'TimeAfterStep_s','ModalFrequency_Hz','LsimFrequency_Hz'}),fullfile(out,'horizon_clarification_trajectory.csv'));
f=fopen(fullfile(out,'horizon_clarification.json'),'w');fprintf(f,'%s\n',jsonencode(result,'PrettyPrint',true));fclose(f);
f=fopen(fullfile(out,'horizon_clarification.txt'),'w');fprintf(f,'POST-RUN HORIZON CLARIFICATION\n\n');
fprintf(f,'%s, case %d remains outside the band at elapsed 29 s: %.12g Hz.\n',result.DesignID,result.CaseID,result.originalFrequencyAtEndpoint_Hz);
fprintf(f,'Continuation to 60 s gives sampled recovery %.6f s and root-refined recovery %.12f s after the step.\n',result.extendedSampledRecovery_s,result.refinedRecovery_s);
fprintf(f,'The two grid-bracketed roots differ by %.4g s. Modal and lsim responses differ by at most %.4g Hz.\n',result.fineGridRootDifference_s,result.modalLsimMaxDifference_Hz);
fprintf(f,'The decreasing modal envelope at 60 s is %.6g Hz, below the 0.001 Hz band. This excludes later re-entry within floating-point arithmetic.\n',tail);
fprintf(f,'Original 30 s metrics remain unchanged: recovery is censored, and the 8 s target fails. No all-target count or eligible matched ranking changes.\n');
fprintf(f,'This is a post-run numerical clarification, not a change to the predeclared protocol or a physical validation.\n');fclose(f);
disp(result);
end
