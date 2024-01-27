using AMDGPU

function mykernel(in)

ix = (workgroupDim().x * (workgroupIdx().x - 1)) + workitemIdx().x
iy = (workgroupDim().y * (workgroupIdx().y - 1)) + workitemIdx().y

if isnan(in[ix,iy])
in[ix,iy] = 0
end

return nothing
end

ac = rand(1024,2048)
ac[3,5] = NaN
ac[1020,2] = NaN
grsize=cld(1024*2048,1024)

ag = ROCArray(ac)
@roc gridsize=(1,2048) groupsize=(1024,1) mykernel(ag)

