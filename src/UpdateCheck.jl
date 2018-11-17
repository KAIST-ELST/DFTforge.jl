
# Check last updated time
import Dates

function getVersionNumber(s)
  version_remote = ""
  for i in 1:length(s)
    v = s[i]
    if occursin("version = ",v)
      a = split(v,[' ','"']; keepempty=false)[3]
      version_remote = VersionNumber(a)
    end
  end
  return version_remote
end

dates_now =  Dates.now()

update_required = false;
#println(pwd())
markerfilePath = joinpath(@__DIR__,"../checkupdate");
versionfilePath = joinpath(@__DIR__,"../Project.toml");
#println(versionfilePath)
if !isfile( markerfilePath)
  update_required = true;
  Base.Filesystem.touch(markerfilePath)
else
  dates_lastchecked = Dates.unix2datetime(Base.Filesystem.mtime("checkupdate"))
  if (dates_now - dates_lastchecked  > Dates.Day(7))
    update_required = true;
  end
end
if update_required
  println(" Update required. Automatic update checking started")
  println(" Visit https://kaist-elst.github.io/DFTforge.jl/ for more update news.")

  try

    # Download metainfo (https://raw.githubusercontent.com/KAIST-ELST/DFTforge.jl/master/Project.toml)
    metainfo_path = "https://github.com/KAIST-ELST/DFTforge.jl/raw/master/Project.toml"
    # If update is required, download and replace
    println(" Checking the updated version ...")
    metadownload_path = download(metainfo_path);

    s_remote = open(metadownload_path) do f
        readlines(f)
    end

    s_local = open(versionfilePath) do f
        readlines(f)
    end



    version_Remote = getVersionNumber(s_remote)
    version_Local = getVersionNumber(s_local)


    if version_Remote > version_Local
      Pkg.update("DFTforge")
      println(" New version: ",version_Remote, " current version: ", version_Local)
      println(" git pull github.com/KAIST-ELST/DFTforge.jl ")
      println(" Or visit https://kaist-elst.github.io/DFTforge.jl/ for update.")
    else
      println(" Current version ",version_Local, " is the newest version.")
    end

  # ziped src https://github.com/KAIST-ELST/DFTforge.jl/archive/master.zip
  catch
    Base.Filesystem.touch(markerfilePath)
  end
end

global DFTforge_VERSION;
function __init__()
  global DFTforge_VERSION;
  s_local = open(versionfilePath) do f
      readlines(f)
  end
  DFTforge_VERSION = getVersionNumber(s_local)
  if 1 == myid()
    println(" DFTforge Version ",string(DFTforge_VERSION))
  end
end

export get_DFTforge_VERSION
function get_DFTforge_VERSION()
  global DFTforge_VERSION;
  return DFTforge_VERSION;
end
