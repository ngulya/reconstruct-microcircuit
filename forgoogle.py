	
gcloud compute scp --zone us-central1-a --recurse  instance-1:~/CloudCircuit/result/ ./

gcloud compute scp --zone us-central1-a  hostfile  instance-1:~/

gcloud compute instances list | awk 'NR>1 {print $8}' > hostfile



gcloud beta compute ssh --zone "us-central1-a" "instance-1" --project "project7-274522"
gcloud beta compute ssh --zone "us-central1-a" "instance-2" --project "project7-274522"



Host 35.193.186.54
  IdentityFile /bckpkeys/google_compute_engine


35.193.110.198 slots=1 max_slots=1
35.202.4.215 slots=1 max_slots=1






gcloud compute --project "project_name" 		disks snapshot "disk_name" --zone "zone_name" --snapshot-names "mysnapshot"
gcloud compute --project steam-kingdom-274515 	disks snapshot instance-1 --zone us-central1-a --snapshot-names mysnapshot

project_name->steam-kingdom-274515
disk_name->instance-1




gcloud compute --project steam-kingdom-274515 disks create instance-2 --zone us-central1-a --source-snapshot mysnapshot --type pd-standard --size 15
gcloud compute --project "project_name" disks create "new_disk_name" --zone "zone_name" --source-snapshot "mysnapshot" --type "pd-standard" --size "10"

new_disk_name->instance-2





gcloud compute --project steam-kingdom-274515 images create instance-3 \\
--source-disk https://www.googleapis.com/compute/v1/projects/steam-kingdom-274515/zones/us-central1-a/disks/instance-2

gcloud compute --project "project_name" images create "new_image" \\
--source-disk https://www.googleapis.com/compute/v1/projects/project_name/zones/zone_name/disks/new_disk_name

new_image->instance-3