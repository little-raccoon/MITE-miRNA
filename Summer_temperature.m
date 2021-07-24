clc
clear all

month = 6;
sum6 = 0;
idealMatrix = int32(10*ones(7680,4717));
for i = 1 :10
    year = 2008+i;
    filename = [num2str(year) '_' num2str(month) '.nc'];
    [LON,LAT,TIME,TMP] = readTMP(filename);
    sum6 = sum6+TMP;
end
cycle6 = round(sum6/-32768);
sum_real6 = (sum6+cycle6*32768)./(idealMatrix-cycle6);

txtValue = sum_real6;
ncid = netcdf.create(['Average_6.nc'],'NC_WRITE');

londim = netcdf.defDim(ncid,'lon',7680);
latdim = netcdf.defDim(ncid,'lat',4717);
timedim = netcdf.defDim(ncid,'time',1);


lonid = netcdf.defVar(ncid,'lon','double',londim);
latid = netcdf.defVar(ncid,'lat','double',latdim);
timeid = netcdf.defVar(ncid,'time','double',timedim);
tmpid = netcdf.defVar(ncid,'tmp','int',[londim,latdim,timedim]);

netcdf.endDef(ncid);


netcdf.putVar(ncid,lonid,LON)
netcdf.putVar(ncid,latid,LAT)
netcdf.putVar(ncid,timeid,6)
netcdf.putVar(ncid,tmpid,txtValue)

netcdf.reDef(ncid);
netcdf.putAtt(ncid,lonid,'long_name','longitude');
netcdf.putAtt(ncid,lonid,'unit','degree');

netcdf.putAtt(ncid,latid,'long_name','latitude');
netcdf.putAtt(ncid,latid,'unit','degree');

netcdf.putAtt(ncid,timeid,'long_name','time');
netcdf.putAtt(ncid,timeid,'unit','2009-2018-06');
netcdf.putAtt(ncid,timeid,'calendar','gregorian');

netcdf.putAtt(ncid,tmpid,'long_name','monthly mean temperature');
netcdf.putAtt(ncid,tmpid,'unit','0.1 degree centigrade');
netcdf.putAtt(ncid,tmpid,'missing_value',-32768);
netcdf.endDef(ncid);
netcdf.close(ncid)
%     dlmwrite([num2str(i) 'yue.txt'],txtValue);


month = 7;
sum7 = 0;
for i = 1 :10
    year = 2008+i;
    filename = [num2str(year) '_' num2str(month) '.nc'];
    [LON1,LAT1,TIME1,TMP1] = readTMP(filename);
    sum7 = sum7+TMP1;
end
cycle7 = round(sum7/-32768);
sum_real7 = (sum7+cycle7*32768)./(idealMatrix-cycle7);


txtValue = sum_real7;
ncid = netcdf.create(['Average_7.nc'],'NC_WRITE');

londim = netcdf.defDim(ncid,'lon',7680);
latdim = netcdf.defDim(ncid,'lat',4717);
timedim = netcdf.defDim(ncid,'time',1);


lonid = netcdf.defVar(ncid,'lon','double',londim);
latid = netcdf.defVar(ncid,'lat','double',latdim);
timeid = netcdf.defVar(ncid,'time','double',timedim);
tmpid = netcdf.defVar(ncid,'tmp','int',[londim,latdim,timedim]);

netcdf.endDef(ncid);


netcdf.putVar(ncid,lonid,LON)
netcdf.putVar(ncid,latid,LAT)
netcdf.putVar(ncid,timeid,7)
netcdf.putVar(ncid,tmpid,txtValue)

netcdf.reDef(ncid);
netcdf.putAtt(ncid,lonid,'long_name','longitude');
netcdf.putAtt(ncid,lonid,'unit','degree');

netcdf.putAtt(ncid,latid,'long_name','latitude');
netcdf.putAtt(ncid,latid,'unit','degree');

netcdf.putAtt(ncid,timeid,'long_name','time');
netcdf.putAtt(ncid,timeid,'unit','2009-2018-07');
netcdf.putAtt(ncid,timeid,'calendar','gregorian');

netcdf.putAtt(ncid,tmpid,'long_name','monthly mean temperature');
netcdf.putAtt(ncid,tmpid,'unit','0.1 degree centigrade');
netcdf.putAtt(ncid,tmpid,'missing_value',-32768);
netcdf.endDef(ncid);
netcdf.close(ncid)
%     dlmwrite([num2str(i) 'yue.txt'],txtValue);


month = 8;
sum8 = 0;
for i = 1 :10
    year = 2008+i;
    filename = [num2str(year) '_' num2str(month) '.nc'];
    [LON2,LAT2,TIME2,TMP2] = readTMP(filename);
    sum8 = sum8+TMP2;
end
cycle8 = round(sum8/-32768);
sum_real8 = (sum8+cycle8*32768)./(idealMatrix-cycle8);


txtValue = sum_real8;
ncid = netcdf.create(['Average_8.nc'],'NC_WRITE');

londim = netcdf.defDim(ncid,'lon',7680);
latdim = netcdf.defDim(ncid,'lat',4717);
timedim = netcdf.defDim(ncid,'time',1);


lonid = netcdf.defVar(ncid,'lon','double',londim);
latid = netcdf.defVar(ncid,'lat','double',latdim);
timeid = netcdf.defVar(ncid,'time','double',timedim);
tmpid = netcdf.defVar(ncid,'tmp','int',[londim,latdim,timedim]);

netcdf.endDef(ncid);


netcdf.putVar(ncid,lonid,LON)
netcdf.putVar(ncid,latid,LAT)
netcdf.putVar(ncid,timeid,8)
netcdf.putVar(ncid,tmpid,txtValue)

netcdf.reDef(ncid);
netcdf.putAtt(ncid,lonid,'long_name','longitude');
netcdf.putAtt(ncid,lonid,'unit','degree');

netcdf.putAtt(ncid,latid,'long_name','latitude');
netcdf.putAtt(ncid,latid,'unit','degree');

netcdf.putAtt(ncid,timeid,'long_name','time');
netcdf.putAtt(ncid,timeid,'unit','2009-2018-08');
netcdf.putAtt(ncid,timeid,'calendar','gregorian');

netcdf.putAtt(ncid,tmpid,'long_name','monthly mean temperature');
netcdf.putAtt(ncid,tmpid,'unit','0.1 degree centigrade');
netcdf.putAtt(ncid,tmpid,'missing_value',-32768);
netcdf.endDef(ncid);
netcdf.close(ncid)
%     dlmwrite([num2str(i) 'yue.txt'],txtValue);

sum678 = 0;
idealMatrix = int32(30*ones(7680,4717));
for i = 1 :10
    for j = 1:3
        month = j+5;
        year = 2008+i;
        filename = [num2str(year) '_' num2str(month) '.nc'];
        [LON123,LAT123,TIME123,TMP123] = readTMP(filename);
        sum678 = sum678+TMP2;
    end
end
cycle678 = round(sum678/-32768);
sum_real678 = (sum678+cycle678*32768)./(idealMatrix-cycle678);


txtValue = sum_real678;
ncid = netcdf.create(['Average_678.nc'],'NC_WRITE');

londim = netcdf.defDim(ncid,'lon',7680);
latdim = netcdf.defDim(ncid,'lat',4717);
timedim = netcdf.defDim(ncid,'time',1);


lonid = netcdf.defVar(ncid,'lon','double',londim);
latid = netcdf.defVar(ncid,'lat','double',latdim);
timeid = netcdf.defVar(ncid,'time','double',timedim);
tmpid = netcdf.defVar(ncid,'tmp','int',[londim,latdim,timedim]);

netcdf.endDef(ncid);


netcdf.putVar(ncid,lonid,LON)
netcdf.putVar(ncid,latid,LAT)
netcdf.putVar(ncid,timeid,678)
netcdf.putVar(ncid,tmpid,txtValue)

netcdf.reDef(ncid);
netcdf.putAtt(ncid,lonid,'long_name','longitude');
netcdf.putAtt(ncid,lonid,'unit','degree');

netcdf.putAtt(ncid,latid,'long_name','latitude');
netcdf.putAtt(ncid,latid,'unit','degree');

netcdf.putAtt(ncid,timeid,'long_name','time');
netcdf.putAtt(ncid,timeid,'unit','2009-2018-678');
netcdf.putAtt(ncid,timeid,'calendar','gregorian');

netcdf.putAtt(ncid,tmpid,'long_name','monthly mean temperature');
netcdf.putAtt(ncid,tmpid,'unit','0.1 degree centigrade');
netcdf.putAtt(ncid,tmpid,'missing_value',-32768);
netcdf.endDef(ncid);
netcdf.close(ncid)
%     dlmwrite([num2str(i) 'yue.txt'],txtValue);



function [LON,LAT,TIME,TMP]=readTMP(filename)

% filename='J:\FTP\20090423test1.hdf';


LON = ncread(filename,'lon');
%  LON=flipud(rot90(LON));

LAT = ncread(filename,'lat');
%  LAT=flipud(rot90(LAT));

TIME = ncread(filename,'time');
%  TIME=flipud(rot90(TIME));

TMP = ncread(filename,'tmp');
%  TMP=flipud(rot90(TMP));
end