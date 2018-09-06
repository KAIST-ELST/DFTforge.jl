
# Check last updated time
import Dates

dates_now =  Dates.now()

update_required = false;

if !isfile("checkupdate")
  update_required = true;
  Base.Filesystem.touch("checkupdate")
else
  dates_lastchecked = Dates.unix2datetime(Base.Filesystem.mtime("checkupdate"))
  if (dates_now - dates_lastchecked  > Dates.Day(7))
    update_required = true;
  end
end
if update_required
  println(" Update required. Automatic update checking started")
  println(" Visit https://kaist-elst.github.io/DFTforge.jl/ for more update news.")
end

try
  # Download metainfo (https://raw.githubusercontent.com/KAIST-ELST/DFTforge.jl/master/Project.toml)
  metainfo_path = "https://github.com/KAIST-ELST/DFTforge.jl/raw/master/Project.toml"
  # If update is required, download and replace
  metadownload_path = download(metainfo_path);

  s_remote = open(metadownload_path) do f
      readlines(f)
  end

  s_local = open("Project.toml") do f
      readlines(f)
  end

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

  version_Remote = getVersionNumber(s_remote)
  version_Local = getVersionNumber(s_local)


  if version_Remote > version_Local
    Pkg.update("DFTforge")
    println(" New version: ",version_Remote, " current version: ", version_Local)
    println(" git pull github.com/KAIST-ELST/DFTforge.jl ")
    println(" Or visit https://kaist-elst.github.io/DFTforge.jl/ for update.")
  end

# ziped src https://github.com/KAIST-ELST/DFTforge.jl/archive/master.zip
finally

end
