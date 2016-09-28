# Loop through all worker nodes, update hosts file and copy ssh public key to it
# The script make the assumption that the node is called %WORKER+<index> and have
# static IP in sequence order

ADMIN_USERNAME=$1
ADMIN_PASSWORD=$2
MASTER_NAME=$3
WORKER_NAME=$4

sudo -u $ADMIN_USERNAME ssh -tt $WORKER_NAME "echo '$ADMIN_PASSWORD' | sudo yum install -y epel-release"
sudo -u $ADMIN_USERNAME ssh -tt $WORKER_NAME "echo '$ADMIN_PASSWORD' | sudo yum install -y sshpass"
sudo -u $ADMIN_USERNAME ssh -tt $WORKER_NAME "echo '$ADMIN_PASSWORD' | mkdir /home/'$ADMIN_USERNAME'/.ssh"
sudo -u $ADMIN_USERNAME ssh -tt $WORKER_NAME "echo '$ADMIN_PASSWORD' | ssh-keygen -f /home/'$ADMIN_USERNAME'/.ssh/id_rsa -t rsa -N ''"
sudo -u $ADMIN_USERNAME ssh -tt $WORKER_NAME "echo '$ADMIN_PASSWORD' | sudo /usr/local/sbin/pbs_mom -p -d /var/spool/torque"

