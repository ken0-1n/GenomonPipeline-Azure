
ADMIN_USERNAME=$1
ADMIN_PASSWORD=$2
MASTER_NAME=$3
WORKER_NAME=$4

sudo -u $ADMIN_USERNAME ssh -tt $WORKER_NAME "echo '$ADMIN_PASSWORD' | sudo /usr/local/sbin/pbs_mom -p -d /var/spool/torque"
