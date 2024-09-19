from celery import shared_task
import os
import docker
from django.conf import settings

CFM_CONFIG_DIR = os.path.join(settings.BASE_DIR, 'VMN', 'config')

@shared_task(bind=True, max_retries=3, default_retry_delay=5)  # 任务绑定到self并允许重试
def run_simulation_task(self, molecule_file_path, user_directory):
    try:
        absolute_user_directory = os.path.abspath(user_directory)
        absolute_molecule_file_path = os.path.abspath(molecule_file_path)

        if not os.path.exists(absolute_molecule_file_path):
            print("Error: molecule.txt not found.")
            return {'status': 'FAILURE', 'error': 'molecule.txt not found.'}

        # 设置Docker容器，限制内存，并设置超时时间
        container = docker.from_env().containers.run(
            image="wishartlab/cfmid",
            command=f"cfm-predict /data/molecule.txt 0.001 /config/param_output.log /config/param_config.txt 0 /data/output.log 0 0",
            volumes={
                absolute_user_directory: {"bind": "/data", "mode": "rw"},
                CFM_CONFIG_DIR: {"bind": "/config", "mode": "ro"}
            },
            platform="linux/amd64",
            detach=True,
            mem_limit="1g",  # 限制内存为1GB
            nano_cpus=1000000000  # 限制CPU为1核，nano_cpus值必须是整数
        )

        # 等待容器任务完成，设置超时时间为600秒
        container.wait(timeout=600)
        logs = container.logs().decode("utf-8")
        print("Container Logs:\n", logs)

        output_file_path = os.path.join(absolute_user_directory, "output.log")
        if os.path.exists(output_file_path):
            return {'status': 'SUCCESS', 'result': output_file_path}  # 返回明确的状态和结果
        else:
            print("output.log file not found.")
            return {'status': 'FAILURE', 'error': 'output.log not found.'}

    except docker.errors.ContainerError as e:
        print(f"Error running container: {e}")
        self.retry(exc=e)  # 使用self.retry来重试

    except Exception as e:
        print(f"An error occurred: {e}")
        self.retry(exc=e)  # 捕获所有其他异常并执行重试

    finally:
        if 'container' in locals():
            container.remove()

    return None





